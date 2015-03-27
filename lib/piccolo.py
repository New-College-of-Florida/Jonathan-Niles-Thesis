"""
Module : Piccolo
Author : Jonathan Niles
Date : October 20, 2014

This module provides a generic loader for the pickle datasets
used to store intermediate information for analysis in the HiC
parsing libraries.
"""

import pickle
import os

PICKLEDIR = "/home/jniles/data/pickles/"

"""Helper functions for reading and writing pickles"""

def _makePath(key):
    """
    """
    path = os.path.join(PICKLEDIR, key + ".pickle")
    return path

def read(path):
    """
    Read a saved pickle file into memory.

    This function takes in an absolute path, verifying
    its existence, and reads it into memory if it exists.
    An exception is thrown if no file exists at that path.

    Parameters
    ----------
    path : String
           Input absolute path
    
    Returns
    -------
    out : Dictionary
          A pickle dictionary with the desired data
    """
    if not os.path.isfile(path):
        raise ValueError("Not a valid path")

    return pickle.load(open(path, 'r'))

def exists(key):
    """
    Determine if a key exists

    This function allows you to check if a pickle exists without raising
    an exception.

    Parameters
    ----------
    key : String
          The key to check existance thereof

    Returns
    -------
    out : Boolean
    """
    return os.path.isfile(_makePath(key))

def load(key):
    """
    Load a previous dataset from a key value

    Parameters
    ---------
    key : String
          The key/id of the data structure

    Returns
    -------
    out : Pickle Dictionary
          The stored data structure or None
    """
    path = _makePath(key)
    data = None

    if os.path.isfile(path):
        data = pickle.load(open(path, 'rb'))
    else:
        print "[WARNING] Path does not exist."
    return data

def save(key, data, overwrite=False):
    """
    Saves a pickable data structure under a key name

    Parameters
    ----------
    key : String
          The key/id for the data structure

    Returns
    -------
    out : None
    """
    path = _makePath(key)
    
    # check if this file exists
    if os.path.isfile(path) and not overwrite:
        raise Exception("File {} already exits")
    
    pickle.dump(data, open(path, 'wb'))

def clean(key):
    """
    Removes a saved pickled data structure

    Paramters
    ---------
    key : String
          The key/id for the data structure

    Returns
    -------
    out : None
    """
    if exists(key):
        path = _makePath(key)
        os.remove(path)
