#returns only function name
def get_function_name():
    l = list()
    with open(r'functions_project.py','r') as fp:
        rld = fp.readlines()
        for i in rld:
            if i[0:3] in 'def':
                try:
                    l.append((i[4:i.index('__')],i[4:-2]))
                except ValueError:
                    l.append((i[4:i.index('(')],i[4:-2]))
    return l

get_function_name()

def get_graph_function_name():
    with open('graph_function.py','r') as fp:
        l = list()
        for i in fp.readlines():
            if i[:3] in 'def':
                print('in')
                l.append((i[4:i.index('(')],i[4:-2]))
    return l

print(get_graph_function_name())

#print(get_function_name())

import functions_project as fp

def run_and_return_function_formula_var(fun):
    vl = []
    rou = '('
    c = 1
    vl = fun[fun.index('(')+1:fun.index(')')].split(',')
    #fun = fun.replace(fun[fun.index('('):],rou)
    for i in range(len(vl)):
        rou += str(i+1)+','
    rou = rou[:-1:]
    rou += ')'
    fun = fun.replace(fun[fun.index('('):fun.index(')')+1],rou)
    #print(fun)
    fv = eval('fp.'+fun)
    try:
        return(fv[2][fv[2].index(':')+1:],tuple(vl))
    except:
        return(fv[2],tuple(vl))

'''
def exp_modiluator(exp,rv):
    lhs,rhs = tuple(exp.split('='))
    p_m_sign_pos = []
    rhs_terms = []
    for i in rhs:
        if i in '+-':
            p_m_sign_pos.append((i,rhs.index(i)))
    print(p_m_sign_pos)
    sp = 0
    a = 0

    for i in range(len(p_m_sign_pos)):
        rhs_terms.append(rhs[sp:p_m_sign_pos[a][1]])
        sp = p_m_sign_pos[a][1] + 1
        a += 1
    rhs_terms.append(rhs[p_m_sign_pos[-1][1]+1:])
    print(rhs_terms)

    for i in p_m_sign_pos:

exp_modiluator('y=mx+cv-u')
'''
