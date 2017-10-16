from chemkin import *

c = chemkin()
c.parse('rxns')

def test_reaction_result():
    assert c.reaction([10,10],[1,2,1],[[1,2,0],[0,0,2]],[[0,0,1],[1,2,0]]) == [-30, -60, 20]

def test_reaction_k_length():
    try:
        c.reaction([10],[1,2,1],[[1,2,0],[0,0,2]],[[0,0,1],[1,2,0]])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_conc_dim():
    try:
        c.reaction([10,10],[1,2],[[1,2,0],[0,0,2]],[[0,0,1],[1,2,0]])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_reac_coef_dim():
    try:
        c.reaction([10,10],[1,2,1],[[1,2],[0,0,2]],[[0,0,1],[1,2,0]])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_prod_coef_dim():
    try:
        c.reaction([10,10],[1,2,1],[[1,2],[0,0,2]],[[0,0],[1,2]])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_prod_coef_complex():
    try:
        c.reaction([10,10],[1,2,1],[[1,2,3],[0,0,2]],[[(0+8j),0,1],[1,2,0]])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_reac_coef_complex():
    try:
        c.reaction([10,10],[1,2,1],[[1,2,3],[0,(1.5+1j),2]],[[1,0,1],[1,2,0]])
    except ValueError as err:
        assert(type(err) == ValueError)



