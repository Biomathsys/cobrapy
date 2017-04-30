# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function

import sympy
import logging

import cobra

logger = logging.getLogger(__name__)


def _is_positive(n):
    """Robustly test if n is positive, yielding True on Exceptions"""
    try:
        if n >= 0:
            return True
        else:
            return False
    except Exception:
        return True


class Frozendict(dict):
    def __init__(self, iterable, **kwargs):
        super(Frozendict, self).__init__(iterable, **kwargs)

    def popitem(self):
        raise AttributeError("'Frozendict' object has no attribute 'popitem")

    def pop(self, k, d=None):
        raise AttributeError("'Frozendict' object has no attribute 'pop")

    def __setitem__(self, key, value):
        raise AttributeError(
            "'Frozendict' object has no attribute '__setitem__")

    def setdefault(self, k, d=None):
        raise AttributeError(
            "'Frozendict' object has no attribute 'setdefault")

    def __delitem__(self, key):
        raise AttributeError(
            "'Frozendict' object has no attribute '__delitem__")

    def __hash__(self):
        return hash(tuple(sorted(self.items())))

    def update(self, E=None, **F):
        raise AttributeError("'Frozendict' object has no attribute 'update")


class AutoVivification(dict):
    """Implementation of perl's autovivification feature. Checkout
    http://stackoverflow.com/a/652284/280182 """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def get_objective_for(model, value):
    """Get or create a reaction for a metabolite or a reaction.

    If value is a Metabolite or a Metabolite id, return any already existing
    demand or exchange reaction.

    Parameters
    ----------
    model : cobra.Model
        The model to for which to get / create a reaction
    value: str, Reaction or Metabolite, None
        A reaction identifier, a Reaction or a Metabolite for which a demand
        reaction is to be created. A None value returns the model's current
        objective.

    Returns
    -------
    model.problem.Objective
    """
    if value is None:
        return model.objective
    if isinstance(value, sympy.Basic):
        return model.problem.Objective(value, direction='max')
    if isinstance(value, model.problem.Objective):
        return value
    try:
        reactions = model.reactions.get_by_any(value)
    except KeyError:
        metabolite = model.metabolites.get_by_any(value)[0]
        reactions = model.reactions.query("^(EX|DM)_{}$".format(metabolite.id))
        if len(reactions) == 0:
            reactions = [model.add_boundary(metabolite, type='demand')]
    return model.problem.Objective(reactions[0].flux_expression)
