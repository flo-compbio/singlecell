"""Functions for working with configuration files for the inDrop pipeline."""

import logging
from pkg_resources import resource_string, resource_filename

import yaml
#from jinja2 import Environment, PackageLoader, select_autoescape

_LOGGER = logging.getLogger(__name__)

#_TEMPLATE_ENV = Environment(
#    loader=PackageLoader('singlecell',
#                         os.path.join('data', 'templates')),
#    autoescape=select_autoescape(['html', 'xml'])
#)


def read_config(config_file):
    """Read the configuration file (in YAML format).

    We rely on a simple configuration file that contains a set of
    sections, with each section containing a list of options. With this simple
    two-level hierarchy, it is straightforward to assign default values to
    options based on the configuration file template included in the package.

    A separate function checks whether mandatory options (e.g., input files)
    are specified.

    TODO: docstring"""

    # read the configuration file template containing the default values for
    # various options
    with open(get_config_template_file(), 'rb') as fh:
        config = yaml.load(fh)

    # read the user-provided configuration file that can be used to override
    # default option values
    with open(config_file, 'rb') as fh:
        user_config = yaml.load(fh)

    errors = False

    for section, entries in user_config.items():
        if section not in config:
            _LOGGER.warning(
                'Ignoring invalid section "%s" in configuration file',
                section)
            continue

        sec = config[section]
        for option, value in entries.items():
            if option not in sec:
                _LOGGER.warning(
                    'Ignoring invalid option "%s" in section "%s" of '
                    'configuration file.', option, section)
            elif value is not None:
                if sec[option] is not None and not \
                        isinstance(value, type(sec[option])):
                    _LOGGER.error(
                        'Wrong data type for option "%s" in section "%s": '
                        'Should be "%s", but got "%s" (%s).',
                        option, section,
                        str(type(sec[option])), str(type(value), str(value)))
                    errors = True
                sec[option] = value

    return config, errors


def write_config(conf, output_file):
    """Write documentation to yaml file."""
    with open(output_file, 'w') as ofh:
        yaml.dump(conf, ofh, default_flow_style=False)


def get_config_template():
    """Create a configuration file template for use with the inDrop pipeline.
    
    TODO: docstring
    """

    return resource_string('singlecell', 'data/indrop/config_template.yaml')\
            .decode('utf-8')


def get_config_template_file():
    """Get the path of the configuration file template.
    
    TODO: docstring
    """
    return resource_filename('singlecell', 'data/indrop/config_template.yaml')
