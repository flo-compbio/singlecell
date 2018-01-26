import pytest

@pytest.fixture(scope='session')
def my_data_pypath(tmpdir_factory):
    """temporary directory for storing test data"""
    pypath = tmpdir_factory.mktemp('singlecell_data', numbered=False)
    return pypath


