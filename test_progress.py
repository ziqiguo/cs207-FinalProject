import progress

def test_progress_result():
    assert progress.progress(10,[1,2,3],[2,1,0]) == 20
    
def test_progress_length():
    try:
        progress.progress(10,[1,2,3],[2,1])
    except ValueError as err:
        assert(type(err) == ValueError)