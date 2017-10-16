from chemkin import *

c = chemkin()

def test_progress_result():
    assert c.progress([10,10],[1,2,1],[[1,2,0],[2,0,2]]) == [40, 10]

def test_progress_k_length():
    try:
        c.progress([10],[1,2,1],[[1,2,0],[2,0,2]])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_coef_dim():
    try:
        c.progress([10,10],[1,2,1],[[1,2],[2,0,2]])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_conc_dim():
    try:
        c.progress([10,10],[1,1],[[1,2,0],[2,0,2]])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_coef_complex():
    try:
        c.progress([10,10],[1,1,1],[[1,(2+1j),0],[2,0,2]])
    except ValueError as err:
        assert(type(err) == ValueError)
