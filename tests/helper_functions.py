import json


def parse_model_parameters(model_parameter_file):
    """
    Parse model parameters from a JSON file.

    Reads the specified JSON file to extract constants, initial conditions,
    and their respective descriptions and units. Returns these values as lists.

    Args:
        model_parameter_file (str): The path to the JSON file containing model parameters.

    Returns:
        constants (list): List of constant values.
        initial_conditions (list): List of initial condition values.
        constant_desc (list): List of constant descriptions with units.
        init_cond_desc (list): List of initial condition descriptions with units.
    """
    constants = []
    initial_conditions = []
    constant_desc = []
    init_cond_desc = []

    with open(model_parameter_file, 'r') as json_file:
        data = json.load(json_file)

    # Populate lists of constants and initial conditions along with their descriptions and units
    for key, value in data['constants'].items():
        constants.append(value['value'])
        constant_desc.append(f"{value['description']} [{value['unit']}]")

    for key, value in data['states'].items():
        initial_conditions.append(value['initial_condition'])
        init_cond_desc.append(f"{value['description']} [{value['unit']}]")

    return constants, initial_conditions, constant_desc, init_cond_desc
