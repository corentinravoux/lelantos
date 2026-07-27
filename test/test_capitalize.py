"""Placeholder / example unit test."""

# test_capitalize.py

def capital_case(x):
    """Return ``x`` with its first character capitalised."""
    return x.capitalize()

def test_capital_case():
    """Check that :func:`capital_case` capitalises the first character."""
    assert capital_case('semaphore') == 'Semaphore'
