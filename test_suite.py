from chemkin import *

def test_progress_system_v2():
    c = chemkin()
    c.parse('rxns.xml')
    c.k_system(1500)
    try:
        c.v2 = [[0., 1., 1., 0., 0., 0+8j],
                [1., 0., 1., 0., 0., 0.],
                [1., 0., 0., 0., 1., 0.]]
        c.progress_system([2., 1., .5, 1., 1., 1.])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_result():
    c = chemkin()
    c.parse('rxns.xml')
    assert c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500) == [-227364086.53073898, 227364586.53073898, 231985198.37073097, -2311055.9199959813, 500.0, -229675142.45073497]

def test_reaction_v1_v2_dim_1():
    c = chemkin()
    c.parse('rxns.xml')
    try:
        c.v1 = [[1,2]]
        c.v2 = [[1,2], [3,4]]
        c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_v1_v2_dim_2():
    c = chemkin()
    c.parse('rxns.xml')
    try:
        c.v1 = [[1,2,3]]
        c.v2 = [[1,2], [3,4]]
        c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)


def test_reaction_prod_coef_complex():
    c = chemkin()
    c.parse('rxns.xml')
    try:
        c.v2 = [[0., 1., 1., 0., 0., 0+8j],
                [1., 0., 1., 0., 0., 0.],
                [1., 0., 0., 0., 1., 0.]]
        c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_concentration_dim():
    c = chemkin()
    c.parse('rxns.xml')
    try:
        c.reaction_rates([2., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        print(err)
        assert(type(err) == ValueError)



print(test_progress_system_v2())
