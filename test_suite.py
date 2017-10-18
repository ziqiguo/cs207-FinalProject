from chemkin import *

fname = 'test_cases/rxns.xml'

def test_k_const():
    c = chemkin()
    try:
        c.k_constant(['a'])
    except TypeError as err:
        assert(type(err) == TypeError)

print(test_k_const())

def test_k_arr_t():
    c = chemkin()
    try:
        c.k_arrhenius(1., 2., -1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_arr_a():
    c = chemkin()
    try:
        c.k_arrhenius(-1., 2., 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_arr_float():
    c = chemkin()
    try:
        c.k_arrhenius(1., 'e', 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_mod_t():
    c = chemkin()
    try:
        c.k_modified(1., 1., 2., -1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_mod_complex():
    c = chemkin()
    try:
        c.k_modified(1., 0+8j, 2., 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_mod_a():
    c = chemkin()
    try:
        c.k_modified(-1., 0., 2., 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_mod_float():
    c = chemkin()
    try:
        c.k_modified(1., 1., 'e', 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_system_k():
    c = chemkin()
    try:
        c.k_system(1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_reac_k0():
    c = chemkin()
    c.parse(fname)
    try:
        c.progress_reaction(0,[1,2,3],[2,1,0])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_reac_result():
    c = chemkin()
    c.parse(fname)
    assert c.progress_reaction(10,[1,2,3],[2,1,0]) == 20

def test_progress_reac_k_length():
    c = chemkin()
    c.parse(fname)
    try:
        c.progress_reaction(10,[1,2,3],[2,1])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_reac_k_zero():
    c = chemkin()
    c.parse(fname)
    try:
        c.progress_reaction(0,[1,2,3],[2,1])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_reac_reac_conc():
    c = chemkin()
    c.parse(fname)
    try:
        c.progress_reaction(10,[0,0,0],[2,1,1])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_reac_reac_coef():
    c = chemkin()
    c.parse(fname)
    try:
        c.progress_reaction(10,[2,1,0],[0,0,0])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_system_k_exist():
    c = chemkin()
    c.parse(fname)
    try:
        c.progress_system([2., 1., .5, 1., 1., 1.])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_system_k_dim():
    c = chemkin()
    c.parse(fname)
    c.k_system(1500)
    try:
        c.k = [1., 2.]
        c.progress_system([2., 1., .5, 1., 1., 1.])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_system_v2():
    c = chemkin()
    c.parse(fname)
    c.k_system(1500)
    try:
        c.v1 = [[0., 1., 1., 0., 0., 0+8j],
                [1., 0., 1., 0., 0., 0.],
                [1., 0., 0., 0., 1., 0.]]
        c.progress_system([2., 1., .5, 1., 1., 1.])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_result():
    c = chemkin()
    c.parse(fname)
    assert c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500) == [-227364086.53073898, 227364586.53073898, 231985198.37073097, -2311055.9199959813, 500.0, -229675142.45073497]

def test_reaction_v1_v2_dim_1():
    c = chemkin()
    c.parse(fname)
    try:
        c.v1 = [[1,2]]
        c.v2 = [[1,2], [3,4]]
        c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_v1_v2_dim_2():
    c = chemkin()
    c.parse(fname)
    try:
        c.v1 = [[1,2,3]]
        c.v2 = [[1,2], [3,4]]
        c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)


def test_reaction_prod_coef_complex():
    c = chemkin()
    c.parse(fname)
    try:
        c.v2 = [[0., 1., 1., 0., 0., 0+8j],
                [1., 0., 1., 0., 0., 0.],
                [1., 0., 0., 0., 1., 0.]]
        c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_concentration_dim():
    c = chemkin()
    c.parse(fname)
    try:
        c.reaction_rates([2., .5, 1., 1., 1.], 1500)
    except ValueError as err:
        print(err)
        assert(type(err) == ValueError)

def test_repr():
    c = chemkin()
    assert(repr(c) == 'chemkin()')
