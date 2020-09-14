import numpy as np


def viterbi(A, O, B, pi, HS, dict):
    """ A - Transition Probabilities
        O - Observations
        B - Emission Probabilities
        pi - Initial Probabilities
        HS - Hidden states (0: Not Gene, 1: Gene)
    """

    n = len(HS)                                   #number of hidden states
    T = O.shape[0]                                #length of observations

    prob_path = np.zeros([n,T], dtype=float)      #tracks path probabilities
    trace_path = np.zeros([n,T], dtype=int)       #tracks best pathing

    for s in range(n):                            #initial states
        prob_path[s,0] = pi[s] * B[s,O[0]]

    for t in range(1,T):                          #iter - each observation
        for s in range(n):                        #iter - each state
            ob_s = B[s,O[t]]                      #P(O[t] | s), emission prob

            options = prob_path[:, t - 1] * A[:, s]     #previous * current-A
            best_s_prev = np.argmax(options)      # best s(t-1) given current s
            best_p_prev = options[best_s_prev] * ob_s   # s(t-1) prob (max)

            prob_path[s,t], trace_path[s,t] = best_p_prev, best_s_prev

    s_final = np.argmax(prob_path[:, -1])         # final state with max prob
    p_final = prob_path[s_final, -1]              # max prob

    path = [s_final]                              # Path - back to front
    s_iter = s_final                              # iter - each states

    for i in range(T - 1, 0, -1):
        s_iter = trace_path[s_iter, i]
        path.append(s_iter)
    path.reverse()                                # Path - front to back

    return path, s_final

    """
    Citation:
    https://web.stanford.edu/~jurafsky/slp3/
    www.cs.jhu.edu/~langmea/resources/lecture_notes/hidden_markov_models.pdf
    """

def main():
    with open("dna_sequence.txt", 'r') as f:
        string = f.read().splitlines()[0]
    genome = list(string)
    dict = {'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3}

    A = np.array([[0.97, 0.03], [0.03, 0.97]])
    O = np.array([dict[i] for i in genome]).reshape(200,1)

    B_0 = np.full((1,4), 0.25)
    B_1 = np.array([0.13, 0.42, 0.22, 0.23]).reshape(1,4)
    B = np.concatenate((B_0, B_1), axis=0) #NG, N @ A, G, C, T respectively

    pi = np.array([0.97, 0.03]).reshape(2,1) #NG vs G

    path, prob = viterbi(A, O, B, pi, [0,1], dict)
    gene_map = np.concatenate((np.array(genome).reshape(1,200), np.array(
        path).reshape(1,200)), axis = 0)

    result = [gene_map[0][i] if gene_map[1][i] == '1' else '-' for i in range(200)]
    result_str = ''.join(result)


    print('Gene composition in nucelotides:', 1 - result.count('-')/200)
    print(result_str)


if __name__ == '__main__':
    main()
