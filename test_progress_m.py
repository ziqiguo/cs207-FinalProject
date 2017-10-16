import progress_m
import progress
def test_progress_m_result():
    assert progress_m.progress_m([10,10],[1,2,1],[[1,2,0],[2,0,2]]) == [40, 10]
    
def test_progress_m_parameters1():
    try:
        progress_m.progress_m([10],[1,2,1],[[1,2,0],[2,0,2]])
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_progress_m_parameters2():
    try:
        progress_m.progress_m([10,10],[1,2,1],[[1,2],[2,0,2]])
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_progress_m_parameters3():
    try:
        progress_m.progress_m([10,10],[1,1],[[1,2,0],[2,0,2]])
    except ValueError as err:
        assert(type(err) == ValueError)
