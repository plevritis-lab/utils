import importlib.util
import sys
from pathlib import Path

import pytest

segmentation_src = Path(__file__).resolve().parent.parent / "1_segmentation" / "src"
celesta_src = Path(__file__).resolve().parent.parent / "2_celesta" / "src"

if str(segmentation_src) not in sys.path:
    sys.path.insert(0, str(segmentation_src))
if str(celesta_src) not in sys.path:
    sys.path.insert(0, str(celesta_src))


def _load_module(name: str, source_directory: Path):
    """loads a module by name from a specific source directory"""
    spec = importlib.util.spec_from_file_location(name, source_directory / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def segmentation_pipeline():
    return _load_module("run_pipeline", segmentation_src)


@pytest.fixture()
def celesta_pipeline():
    return _load_module("run_pipeline", celesta_src)
