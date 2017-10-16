from chemkin import *

c = chemkin()

def test_progress_u_result():
    assert c.progress_u(10,[1,2,3],[2,1,0]) == 20

def test_progress_u_k_length():
    try:
        c.progress_u(10,[1,2,3],[2,1])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_u_k_zero():
    try:
        c.progress_u(0,[1,2,3],[2,1])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_u_reac_conc():
    try:
        c.progress_u(10,[0,0,0],[2,1,1])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_u_reac_coef():
    try:
        c.progress_u(10,[2,1,0],[0,0,0])
    except ValueError as err:
        assert(type(err) == ValueError)
