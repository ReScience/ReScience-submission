# -*- coding: utf-8 -*-

from configparser import ConfigParser


def get_units_section(config):
    """
    Get config file options from section containing strings.
    
    :param config: ConfigParser object.
    :param section: Name of the section to be read.
    """
    section = 'Units'
    options = config.options(section)
    section_dict = {}    
    for option in options:
        value = config.get(section, option).split('*')
        if len(value) == 2:
            value = float(value[0]) * float(section_dict[value[1].strip()])
        elif len(value) == 4:
            unit = value[1].split('/')
            unit = float(section_dict[unit[0].strip()]) / float(section_dict[unit[1].strip()])**float(value[-1])
            value = float(value[0]) * unit
        else:
            value = value[0]
            try:
                value = float(value)
            except ValueError:
                value = value.split('/')
                value = float(section_dict[value[0].strip()]) * float(section_dict[value[1].strip()])
        section_dict[option] = value
    return section_dict


def get_param_section(config, units):
    """
    Get config file options from section containing strings.
    
    :param config: ConfigParser object.
    :param section: Name of the section to be read.
    """
    section = 'Parameter'
    options = config.options(section)
    section_dict = {}    
    for option in options:
        value = config.get(section, option).split('*')
        if len(value) == 1:
            value = float(value[0])
        elif len(value) == 2:
            unit = value[1].split('/')
            try:
                unit = units[unit[0].strip()] / units[unit[1].strip()]
            except IndexError:
                unit = units[unit[0].strip()]
            except KeyError:
                unit = 1 / units[unit[1].strip()]
            value = float(value[0]) * unit
        elif len(value) == 3:
            subunit = value[1].split('/')
            try:
                unit = units[subunit[0].strip()] / units[subunit[1].strip()] / units[value[2].strip()]
            except KeyError:
                unit = 1 / units[subunit[1].strip()] / units[value[2].strip()]
            value = float(value[0]) * unit
        elif len(value) == 4:
            unit = value[1].split('/')
            try:
                unit = units[unit[0].strip()] / units[unit[1].strip()]**int(value[3])
            except IndexError:
                unit = units[unit[0].strip()]**int(value[3])
            value = float(value[0]) * unit
        section_dict[option] = value
    return section_dict
    

def read_config(fname):
    """
    Reads config.cfg file.
        
    Reads configuration file and sets up parameters for the simulation.
    
    :param fname: Filename of the configuration file.
    """
    config = ConfigParser()
    config.optionxform = str 
    config.read(fname)
    # Files
    unts = get_units_section(config)
    # Arteries
    param = get_param_section(config, unts)
    return unts, param


def read_units(fname):
    """
    Reads config.cfg file.
        
    Reads configuration file and sets up parameters for the simulation.
    
    :param fname: Filename of the configuration file.
    """
    config = ConfigParser()
    config.optionxform = str 
    config.read(fname)
    # Files
    unts = get_units_section(config)
    return unts