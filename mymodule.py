# mymodule.py

"""
This module provides basic mathematical and greeting functions.

It includes:
- add: a function to add two numbers.
- greet: a function to return a greeting message.
"""

def add(a, b):
    """
    Return the sum of two numbers.
    
    Parameters:
    a (int, float): The first number.
    b (int, float): The second number.

    Returns:
    int, float: The sum of a and b.
    """
    return a + b

def greet(name):
    """
    Return a greeting string for the given name.
    
    Parameters:
    name (str): The name to greet.

    Returns:
    str: A greeting message.
    """
    return f"Hello, {name}!"
