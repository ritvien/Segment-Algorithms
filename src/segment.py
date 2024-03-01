
def slidingwindowsegment(sequence, create_segment, compute_error, max_error, seq_range=None):
    """
    Return a list of line segments that approximate the sequence.

    The list is computed using the sliding window technique. 

    Parameters
    ----------
    sequence : sequence to segment
    create_segment : a function of two arguments (sequence, sequence range) that returns a line segment that approximates the sequence data in the specified range
    compute_error: a function of two argments (sequence, segment) that returns the error from fitting the specified line segment to the sequence data
    max_error: the maximum allowable line segment fitting error

    """
    if not seq_range:
        seq_range = (0,len(sequence)-1)

    start = seq_range[0]
    end = start
    result_segment = create_segment(sequence,(seq_range[0],seq_range[1]))
    while end < seq_range[1]:
        end += 1
        test_segment = create_segment(sequence,(start,end))
        error = compute_error(sequence,test_segment)
        if error <= max_error:
            result_segment = test_segment
        else:
            break

    if end == seq_range[1]:
        return [result_segment]
    else:
        return [result_segment] + slidingwindowsegment(sequence, create_segment, compute_error, max_error, (end-1,seq_range[1]))
        
def bottomupsegment(sequence, create_segment, compute_error, max_error):
    """
    Return a list of line segments that approximate the sequence.
    
    The list is computed using the bottom-up technique.
    
    Parameters
    ----------
    sequence : sequence to segment
    create_segment : a function of two arguments (sequence, sequence range) that returns a line segment that approximates the sequence data in the specified range
    compute_error: a function of two argments (sequence, segment) that returns the error from fitting the specified line segment to the sequence data
    max_error: the maximum allowable line segment fitting error
    
    """
    segments = [create_segment(sequence,seq_range) for seq_range in zip(range(len(sequence))[:-1],range(len(sequence))[1:])]
    mergesegments = [create_segment(sequence,(seg1[0],seg2[2])) for seg1,seg2 in zip(segments[:-1],segments[1:])]
    mergecosts = [compute_error(sequence,segment) for segment in mergesegments]
    while min(mergecosts) < max_error:
            idx = mergecosts.index(min(mergecosts))
            segments[idx] = mergesegments[idx]
            del segments[idx+1]

            if idx > 0:
                mergesegments[idx-1] = create_segment(sequence,(segments[idx-1][0],segments[idx][2]))
                mergecosts[idx-1] = compute_error(sequence,mergesegments[idx-1])

            if idx+1 < len(mergecosts):
                mergesegments[idx+1] = create_segment(sequence,(segments[idx][0],segments[idx+1][2]))
                mergecosts[idx+1] = compute_error(sequence,mergesegments[idx])

            del mergesegments[idx]
            del mergecosts[idx]
            if mergecosts==[]:
                break                
    return segments
    
    
def topdownsegment(sequence, create_segment, compute_error, max_error, seq_range=None):
    """
    Return a list of line segments that approximate the sequence.
    
    The list is computed using the bottom-up technique.
    
    Parameters
    ----------
    sequence : sequence to segment
    create_segment : a function of two arguments (sequence, sequence range) that returns a line segment that approximates the sequence data in the specified range
    compute_error: a function of two argments (sequence, segment) that returns the error from fitting the specified line segment to the sequence data
    max_error: the maximum allowable line segment fitting error
    
    """
    if not seq_range:
        seq_range = (0,len(sequence)-1)

    bestlefterror,bestleftsegment = float('inf'), None
    bestrighterror,bestrightsegment = float('inf'), None
    bestidx = None

    for idx in range(seq_range[0]+1,seq_range[1]):
        segment_left = create_segment(sequence,(seq_range[0],idx))
        error_left = compute_error(sequence,segment_left)
        segment_right = create_segment(sequence,(idx,seq_range[1]))
        error_right = compute_error(sequence, segment_right)
        if error_left + error_right < bestlefterror + bestrighterror:
            bestlefterror, bestrighterror = error_left, error_right
            bestleftsegment, bestrightsegment = segment_left, segment_right
            bestidx = idx
    
    if bestlefterror <= max_error:
        leftsegs = [bestleftsegment]
    else:
        leftsegs = topdownsegment(sequence, create_segment, compute_error, max_error, (seq_range[0],bestidx))
    
    if bestrighterror <= max_error:
        rightsegs = [bestrightsegment]
    else:
        rightsegs = topdownsegment(sequence, create_segment, compute_error, max_error, (bestidx,seq_range[1]))
    
    return leftsegs + rightsegs

def SWABsegment(sequence, create_segment, compute_error, max_error, seg_num=5):
    """
    Return a list of line segments that approximate the sequence using the SWAB algorithm.

    Parameters
    ----------
    max_error : maximum allowable line segment fitting error
    seg_num : integer indicating the number of segments (5 or 6)
    sequence : input sequence of data points
    sequence_right : sequence that hasn't been segment, the window will slide to.
    Returns
    -------
    Seg_TS : list of line segments approximating the sequence
    """
    Seg_TS = []  

    w = int(len(sequence)/seg_num)
    lower_bound = w // 2
    upper_bound = 2 * w
    start = 0
    sequence_right = sequence 
    w_in = 0
    while sequence_right:
        T = bottomupsegment(sequence[start:start+w], create_segment, compute_error, max_error)
        w_out = T[0][2]
        w = w-w_out   #w = w -w'

        (x0,y0,x1,y1) = create_segment(sequence,(start,start+w_out))
        start += w_out
        Seg_TS.append((x0,y0,x1,y1))
        
        sequence_right = sequence[start:]
        
        w_in = len(BEST_LINE(max_error, sequence_right,compute_error))
        w += w_in
        if w < lower_bound:
            w = upper_bound
        elif w > upper_bound:
            w = lower_bound
        if start == len(sequence)-1:
            return Seg_TS
    
    if start == len(sequence)-1:
        return Seg_TS
    else:
        print('Last bottom-up:',len(sequence_right))
        last_bottom = bottomupsegment(sequence_right, create_segment, compute_error, max_error)
        last_add = []
        for seg in last_bottom:
            x0 = seg[0] + start
            y0 = seg[1]
            x1 = seg[2] + start
            last_add.append((x0,y0,x1,y1))
        return Seg_TS + last_add


def BEST_LINE(max_error, sequence,compute_error):
    """
    Return a segment S to approximate the next potential segment with error ≤ max_error.

    Parameters
    ----------
    max_error : maximum allowable line segment fitting error
    sequence : input sequence of data points

    Returns
    -------
    S : segment to approximate the next potential segment with error ≤ max_error
    """
    S = [0,sequence[0]]
    for i in range(len(sequence)-1):
        S.append(i+1)
        S.append(sequence[i+1])
        error = compute_error(sequence, (0,sequence[0],i+1,sequence[i+1]))
        if error > max_error:
            S.pop()
            break
    
    return S


def get_ts_segments(segment_method,sequence, create_segment, compute_error, max_error):
    """
    Return a list of line segments that approximate the sequence
    and return a list of time series segments
    
    Parameters
    ----------
    segment_method: sliding windows, topdown or buttom up
    sequence : sequence to segment
    create_segment : a function of two arguments (sequence, sequence range) that returns a line segment that approximates the sequence data in the specified range
    compute_error: a function of two argments (sequence, segment) that returns the error from fitting the specified line segment to the sequence data
    max_error: the maximum allowable line segment fitting error
    
    """
    approximate_segments = segment_method(sequence, create_segment, compute_error, max_error)
    # draw_segments(segments)
    ts_segments = []
    for seg in approximate_segments:
        x0 = seg[0]
        x1 = seg[2]
        ts_seg = sequence[x0:x1]
        ts_segments.append(ts_seg)
        
    return ts_segments, approximate_segments
