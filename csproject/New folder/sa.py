def add_fun():
    with open('functions_project.py','a') as fp:
        fn = input('function name : ').replace(' ','_')
        ln = input('function last name : ')
        fo = input('formula : ')
        lhs,rhs = tuple(fo.split('='))
        v = ''
        for i in rhs:
            if i.isalpha() and i not in v:
                v = v + str(i) +','
        v = v[:-1]
        fun = f"\ndef {fn}__{ln}({v}):\n\treturn ('{lhs}=',{rhs},':{fo}')\n"
        fp.write(fun)
    print(fun)
    print('\n','-'*5,'new','-'*5,'\n')
