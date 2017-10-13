import reaction
def test_reaction_result():
    assert reaction.reaction([10,10],[1,2,1],[[1,2,0],[0,0,2]],[[0,0,1],[1,2,0]]) == [-30, -60, 20]
    
def test_reaction_parameters1():
    try:
        reaction.reaction([10],[1,2,1],[[1,2,0],[0,0,2]],[[0,0,1],[1,2,0]])
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_reaction_parameters2():
    try:
        reaction.reaction([10,10],[1,2],[[1,2,0],[0,0,2]],[[0,0,1],[1,2,0]])
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_reaction_parameters3():
    try:
        reaction.reaction([10,10],[1,2,1],[[1,2],[0,0,2]],[[0,0,1],[1,2,0]])
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_reaction_parameters4():
    try:
        reaction.reaction([10,10],[1,2,1],[[1,2,0],[0,0,2]],[[0,0],[1,2]])
    except ValueError as err:
        assert(type(err) == ValueError)