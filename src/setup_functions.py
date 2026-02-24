# ---------------------------------
# wisdm - environment setup
# ApexRMS, March 2024
# ---------------------------------

# Shared conda/GDAL/PROJ environment setup for all wisdm Python scripts.
# This module uses only stdlib so it can be safely imported before any
# non-stdlib packages have been loaded.

import glob
import os
import platform
import sys

# ---- Temporary debug flag -- set to False to silence debug output ----------
_DEBUG = True

# Buffer for messages that arrive before pysyncrosim is importable.
_debug_messages = []


def _dbg(msg):
    if not _DEBUG:
        return
    formatted = f"[setup_functions] {msg}"
    ps = sys.modules.get("pysyncrosim")
    if ps is not None:
        try:
            ps.environment.update_run_log(formatted)
            return
        except Exception:
            pass
    _debug_messages.append(formatted)


# Execute immediately on import: remove user site-packages before any
# non-stdlib packages are imported. This is critical on Windows where a
# user-installed pysyncrosim or rasterio (e.g. for a different Python version)
# can shadow the conda env packages and cause DLL version conflicts.
_before = [p for p in sys.path if "AppData\\Roaming\\Python" in p
           or "AppData/Roaming/Python" in p]
sys.path[:] = [p for p in sys.path if not (
    "AppData\\Roaming\\Python" in p or "AppData/Roaming/Python" in p)]
if _before:
    _dbg(f"Removed user site-packages from sys.path: {_before}")
else:
    _dbg("No user site-packages found in sys.path")


def setup_conda_env():
    """Ensure the active conda environment's packages take priority over system
    Python packages. Must be called before importing any non-stdlib packages.

    Handles three cases:
    - Conda already activated via 'conda activate' or 'conda run' (CONDA_SHLVL >= 1):
      skips PATH manipulation since conda already set everything up.
    - CONDA_PREFIX is set but conda was not activated through the shell:
      manually mirrors what 'conda activate' does to PATH and DLL directories.
    - CONDA_PREFIX is not set: attempts to derive the conda prefix from
      sys.executable (e.g. when Python is invoked directly by SyncroSim).
    """
    _dbg(f"CONDA_PREFIX: {os.environ.get('CONDA_PREFIX')}")
    _dbg(f"CONDA_SHLVL:  {os.environ.get('CONDA_SHLVL', '(not set)')}")
    _dbg(f"sys.executable: {sys.executable}")

    # If conda already activated the environment, its PATH is correct; leave it.
    conda_shlvl = os.environ.get("CONDA_SHLVL", "0")
    if conda_shlvl != "0":
        _dbg(f"CONDA_SHLVL={conda_shlvl}: conda already activated -- "
             "skipping PATH and sys.path manipulation")
        _log_gdal_on_path()
        return

    # Detect conda prefix from env var, or fall back to deriving from sys.executable
    # (happens when Python is called directly without conda activation).
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        python_path = os.path.dirname(sys.executable)
        _dbg(
            f"CONDA_PREFIX not set, checking sys.executable path: {python_path}")
        if "envs" in python_path or "conda" in python_path.lower():
            conda_prefix = python_path
            _dbg(f"Derived conda_prefix from sys.executable: {conda_prefix}")
        else:
            _dbg("Could not derive conda_prefix from sys.executable")

    if not conda_prefix or not os.path.exists(conda_prefix):
        _dbg(
            f"conda_prefix not found or does not exist: {conda_prefix!r} -- skipping setup")
        return

    # Ensure CONDA_PREFIX is set in the environment (some packages expect it).
    if "CONDA_PREFIX" not in os.environ:
        os.environ["CONDA_PREFIX"] = conda_prefix
        _dbg(f"Set CONDA_PREFIX env var to: {conda_prefix}")

    _dbg("Running manual conda path setup (conda not activated via shell)")

    if platform.system() == "Windows":
        paths_to_add = [
            os.path.join(conda_prefix, "Library", "mingw-w64", "bin"),
            os.path.join(conda_prefix, "Library", "usr", "bin"),
            os.path.join(conda_prefix, "Library", "bin"),
            os.path.join(conda_prefix, "Scripts"),
            os.path.join(conda_prefix, "bin"),
            conda_prefix,
        ]
        for path in paths_to_add:
            if os.path.exists(path):
                os.environ["PATH"] = path + os.pathsep + \
                    os.environ.get("PATH", "")
                if hasattr(os, "add_dll_directory"):
                    try:
                        os.add_dll_directory(path)
                        _dbg(f"  add_dll_directory: {path}")
                    except (FileNotFoundError, OSError) as e:
                        _dbg(f"  add_dll_directory FAILED for {path}: {e}")
            else:
                _dbg(f"  path does not exist, skipping: {path}")

    # Add the conda environment's site-packages to sys.path.
    # Try Windows layout (Lib) first, then Unix layout (lib/pythonX.Y).
    conda_site_packages = os.path.join(conda_prefix, "Lib", "site-packages")
    if not os.path.exists(conda_site_packages):
        python_version = f"python{sys.version_info.major}.{sys.version_info.minor}"
        conda_site_packages = os.path.join(
            conda_prefix, "lib", python_version, "site-packages")

    if os.path.exists(conda_site_packages) and conda_site_packages not in sys.path:
        sys.path.insert(0, conda_site_packages)
        _dbg(
            f"Inserted conda site-packages into sys.path: {conda_site_packages}")
    else:
        _dbg(
            f"conda site-packages already in sys.path or not found: {conda_site_packages}")

    _log_gdal_on_path()


def _log_gdal_on_path():
    """Log which directories on PATH contain GDAL DLLs."""
    if not _DEBUG or platform.system() != "Windows":
        return
    gdal_dirs = []
    for p in os.environ.get("PATH", "").split(os.pathsep):
        if p and os.path.isdir(p) and glob.glob(os.path.join(p, "gdal*.dll")):
            dlls = [os.path.basename(f)
                    for f in glob.glob(os.path.join(p, "gdal*.dll"))]
            gdal_dirs.append((p, dlls))
    if gdal_dirs:
        _dbg(f"GDAL DLLs found in PATH ({len(gdal_dirs)} location(s)):")
        for path, dlls in gdal_dirs:
            _dbg(f"  {path}: {dlls}")
    else:
        _dbg("No GDAL DLLs found in PATH")


def check_gdal_version():
    """On Windows, when multiple GDAL DLL installations are detected on PATH,
    remove any paths that contain a GDAL version older than 3.6 to prevent
    conflicts. No-op on non-Windows platforms or if win32api is unavailable.
    """
    if platform.system() != "Windows":
        return

    try:
        from win32api import GetFileVersionInfo, HIWORD, LOWORD

        gdal_installations = []
        for p in os.environ.get("PATH", "").split(os.pathsep):
            if p and glob.glob(os.path.join(p, "gdal*.dll")):
                gdal_installations.append(os.path.abspath(p))

        _dbg(
            f"check_gdal_version: {len(gdal_installations)} GDAL installation(s) on PATH")

        if len(gdal_installations) > 1:
            for folder in gdal_installations:
                filenames = [f for f in os.listdir(folder)
                             if f.startswith("gdal") and f.endswith(".dll")]
                for filename in filenames:
                    filepath = os.path.join(folder, filename)
                    if not os.path.exists(filepath):
                        os.environ["PATH"] = os.pathsep.join(
                            [p for p in os.environ["PATH"].split(os.pathsep)
                             if folder not in p])
                        continue
                    try:
                        info = GetFileVersionInfo(filepath, "\\")
                    except Exception:
                        continue
                    major_version = HIWORD(info["FileVersionMS"])
                    minor_version = LOWORD(info["FileVersionMS"])
                    _dbg(f"  {filename}: version {major_version}.{minor_version} "
                         f"in {folder}")
                    if major_version < 3 or minor_version < 6:
                        _dbg(
                            f"  -> removing outdated GDAL path from PATH: {folder}")
                        os.environ["PATH"] = os.pathsep.join(
                            [p for p in os.environ["PATH"].split(os.pathsep)
                             if folder not in p])
        elif len(gdal_installations) == 1:
            _dbg(f"  Single GDAL installation, no conflict check needed: "
                 f"{gdal_installations[0]}")

    except ImportError:
        _dbg("win32api not available; skipping GDAL version check")


def setup_gdal_proj(myLibrary):
    """Set GDAL and PROJ environment variables for the active conda environment.
    Must be called after connecting to SyncroSim (myLibrary must be available).

    Returns a dict of environment variables to forward to Dask workers, or
    None if conda is not in use.
    """
    import pysyncrosim as ps
    import pyproj

    # Flush messages buffered before pysyncrosim was available.
    global _debug_messages
    for msg in _debug_messages:
        try:
            ps.environment.update_run_log(msg)
        except Exception:
            pass
    _debug_messages = []

    if myLibrary.datasheets("core_Option").UseConda.item() != "Yes":
        _dbg("UseConda=No -- skipping GDAL/PROJ env var setup")
        return None

    conda_env_path = os.environ.get("CONDA_PREFIX") or sys.prefix
    _dbg(f"setup_gdal_proj: conda_env_path={conda_env_path}")

    if platform.system() == "Windows":
        library_folder = os.path.join(conda_env_path, "Library")
    else:
        library_folder = conda_env_path

    gdal_folder = os.path.join(library_folder, "share", "gdal")
    proj_folder = os.path.join(library_folder, "share", "proj")
    certifi_folder = os.path.join(library_folder, "ssl", "cacert.pem")

    ps.environment.update_run_log("GDAL path: " + gdal_folder)
    ps.environment.update_run_log("PROJ path: " + proj_folder)
    os.environ["GDAL_DATA"] = gdal_folder
    os.environ["GDAL_CURL_CA_BUNDLE"] = certifi_folder
    os.environ["PROJ_DATA"] = proj_folder
    os.environ["PROJ_CURL_CA_BUNDLE"] = certifi_folder
    os.environ["PROJ_LIB"] = proj_folder
    pyproj.datadir.set_data_dir(proj_folder)
    pyproj.network.set_ca_bundle_path(certifi_folder)
    ps.environment.update_run_log(
        "pyproj data directory: " + pyproj.datadir.get_data_dir())

    # Return env vars for Dask workers to inherit.
    worker_env = {"GDAL_DATA": gdal_folder, "PROJ_LIB": proj_folder}
    if os.path.exists(certifi_folder):
        worker_env["GDAL_CURL_CA_BUNDLE"] = certifi_folder
        worker_env["PROJ_CURL_CA_BUNDLE"] = certifi_folder

    return worker_env
