# ---------------------------------
# wisdm - Python environment smoke tests
# ApexRMS, October 2026
# ---------------------------------

# Checks that the wisdm Python conda environment matches src/wisdm-py-conda.yml
# (every pinned package imports, at the pinned version) and that the spatial
# libraries GDAL and PROJ work. Run these after changing a version in the yml
# file. They are skipped outside the wisdm Python conda environment.
#
# The spatial packages are imported inside each test rather than at the top of
# the file, so that running the tests outside the conda environment skips them
# instead of failing on a missing package.

import importlib
import importlib.metadata
import os
import platform

import numpy as np
import pytest

import setup_functions as sf

condaFile = os.path.join(
    os.path.dirname(__file__), "..", "..", "src", "wisdm-py-conda.yml")


def readCondaFile(path):
    """Return the environment name and a {package: version} dict of pins.

    Lines like "  - rasterio=1.3.10" become {"rasterio": "1.3.10"}.
    """
    envName = None
    pinnedVersions = {}
    with open(path, encoding="utf-8") as condaYml:
        for line in condaYml:
            line = line.strip()
            if line.startswith("name:"):
                envName = line.split(":", 1)[1].strip()
            elif line.startswith("- ") and "=" in line:
                packageName, version = line[2:].split("=", 1)
                pinnedVersions[packageName.strip()] = version.strip()
    return envName, pinnedVersions


envName, pinnedVersions = readCondaFile(condaFile)
activeEnvName = os.path.basename(
    os.path.normpath(os.environ.get("CONDA_PREFIX", "")))

pytestmark = pytest.mark.skipif(
    activeEnvName != envName,
    reason=f"only runs inside the {envName} conda environment")


@pytest.fixture(scope="module", autouse=True)
def condaPaths():
    """Set up conda paths the same way the wisdm transformers do before
    importing spatial packages."""
    sf.setupCondaEnv()


def test_python_is_pinned_version():
    assert platform.python_version() == pinnedVersions["python"]


# Each pinned package's import name matches its conda name. If a future pin
# doesn't (e.g. a package installed as "foo-bar" but imported as "foo"), give
# it its own test.
@pytest.mark.parametrize(
    "packageName", [name for name in pinnedVersions if name != "python"])
def test_pinned_package_imports_at_pinned_version(packageName):
    importlib.import_module(packageName)

    assert importlib.metadata.version(packageName) == \
        pinnedVersions[packageName]


def test_rasterio_and_rioxarray_read_a_geotiff(tmp_path):
    """Checks that GDAL works."""
    import rasterio
    import rasterio.transform
    import rioxarray

    rasterPath = tmp_path / "smoke-test.tif"
    values = np.arange(100, dtype="float32").reshape(10, 10)
    profile = {
        "driver": "GTiff",
        "height": 10,
        "width": 10,
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": rasterio.transform.from_origin(0, 10, 1, 1),
    }

    with rasterio.open(rasterPath, "w", **profile) as rasterFile:
        rasterFile.write(values, 1)

    with rasterio.open(rasterPath) as rasterFile:
        assert rasterFile.read(1).sum() == values.sum()
        assert rasterFile.crs.to_epsg() == 4326

    raster = rioxarray.open_rasterio(rasterPath)
    assert raster.shape == (1, 10, 10)
    raster.close()

    assert rasterio.__gdal_version__


def test_geopandas_reprojects_a_point():
    """Checks that PROJ works."""
    import geopandas as gpd

    points = gpd.GeoSeries(gpd.points_from_xy([-120], [40]), crs="EPSG:4326")

    projected = points.to_crs("EPSG:3857")

    # Web Mercator x = longitude * 20037508.34 / 180 = -13358338.9
    assert projected.x[0] == pytest.approx(-13358338.9, rel=1e-6)


def test_dask_sums_an_array():
    import dask.array as da

    assert da.ones((10, 10), chunks=5).sum().compute() == 100
