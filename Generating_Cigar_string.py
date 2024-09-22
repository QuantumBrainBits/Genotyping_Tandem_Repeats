# Function to generate the Cigar string from given, refn seqn and read seqn. (Global alignment)

def Generate_cigar(reference_seq, query_seq):

    # Extract CIGAR string from the first alignment
    query_lindex  = len(query_seq) - len(query_seq.lstrip('-')) # query indexs helps us to remove "leading spaces" & "trailing spaces"
    query_rindex  = len(query_seq) - len(query_seq.rstrip('-')) 
    
    
    # appending the cigar ops in order.
    cigar_ops = []
    
    for pos in range(0+query_lindex, len(query_seq)-query_rindex):
    
        # refn and query seq each base to check its type.
        rbase = reference_seq[pos]
        qbase = query_seq[pos]
        
        # checking "Macthes": represents "M"
        if rbase == qbase :
            cigar_ops.append('M')
    
        # checking "Insertion": "I"
        elif rbase == '-':
            cigar_ops.append('I')
            
        # checking "Deletion": "D"
        elif qbase == '-':
            cigar_ops.append('D')
        
        # checking "Sub": "X"
        elif rbase != qbase:
            cigar_ops.append('X')
            
    
    # list cigar ops into joint cigar string.
    cigar_ops_values = {'M':0, 'I':0, 'D':0, 'X':0}
    cigar_type     = []
    cigar_num      = []
    
    current_ops = cigar_ops[0]
    for i,ops in enumerate(cigar_ops):
    
        # keeps adding the type ops onl if consicutive.
        if ops == current_ops:
            cigar_ops_values[ops] += 1
        
        else: # here we realised order is changed with new ops.
            cigar_type.append(current_ops) # append the current info
            cigar_num.append(cigar_ops_values[current_ops])
            cigar_ops_values.update({}.fromkeys(cigar_ops_values,0)) # its necessary to keep the track of cigar ops so we clear dict so we only update order wise.
            cigar_ops_values[ops] += 1
            current_ops = ops
            
        if ops == current_ops and len(cigar_ops)-1 == int(i):
            cigar_type.append(current_ops)
            cigar_num.append(cigar_ops_values[current_ops])
        

    # return the cigar_type and cigar_num bases.
    query_seq1 = query_seq[0+query_lindex : len(query_seq)-query_rindex]
    return [cigar_type, cigar_num, query_seq1]