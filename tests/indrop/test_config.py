"""Test functions for the `indrop.config` module."""

import logging
from pkg_resources import resource_filename

# import pytest

from singlecell.indrop import config

_LOGGER = logging.getLogger(__name__)


def test_get_config_template():
    """Tests getting the configuration template."""

    conf_str = config.get_config_template()
    assert isinstance(conf_str, str)


def test_read_config():
    """Tests reading the configuration from file."""

    config_file = resource_filename(
        'singlecell', 'data/indrop/config_template.yaml')

    conf = config.read_config(config_file)
    assert isinstance(conf, dict)
