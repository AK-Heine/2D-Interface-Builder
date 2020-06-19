from interfacebuilder.main import interface

import time
import pytest


def test_heterostack():
    A = "tests/graphene.xyz"
    B = "tests/MoS2_2H_1L.xyz"
    if1 = interface(
        A,
        B,
        N_translations=10,
        angle_limits=(0, 30),
        angle_stepsize=0.5,
        tolerance=0.05,
    )
    assert len(if1.pairs) > 0, "Pairs not found"
    if1.analyze_results(weight=0.5, distance=3, optimize=False)
    assert len(if1.solved) > 0, "Structures not build"
    if1.analyze_results(weight=0.5, distance=3, optimize=True)
    assert len(if1.solved) > 0, "Opt not working"
