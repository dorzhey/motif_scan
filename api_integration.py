import requests

def query_jaspar_api(species, pwm=None):
    """
    Queries the JASPAR API for motifs.

    Args:
        species (str): Species name (e.g., "human").
        pwm (str, optional): Position Weight Matrix identifier.

    Returns:
        list: A list of motifs retrieved from the database.
    """
    base_url = "https://jaspar.genereg.net/api/v1/matrix/"
    params = {'species': species}
    if pwm:
        params['pwm'] = pwm

    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Failed to query JASPAR API: {response.status_code} {response.text}")

def query_hocomoco_api(species, pwm=None):
    """
    Queries the HOCOMOCO API for motifs.

    Args:
        species (str): Species name (e.g., "human").
        pwm (str, optional): Position Weight Matrix identifier.

    Returns:
        list: A list of motifs retrieved from the database.
    """
    base_url = "https://hocomoco.autosome.org/api/v1/matrix/"
    params = {'species': species}
    if pwm:
        params['pwm'] = pwm

    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Failed to query HOCOMOCO API: {response.status_code} {response.text}")
