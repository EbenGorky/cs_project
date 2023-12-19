from tkinter import *
from tkinter import messagebox
import oxilanari_functions as oxi
import functions_project as f_p
import graph_function as g_p 

oxi_info = oxi.get_function_name()
oxi_gra_info = oxi.get_graph_function_name()
#print(oxi_info)

global win
win = Tk()
win.title('CS project')
win.geometry('600x600')

fram_intro = None
fram_for = None
fram_gra_con = None
fram_ = None
fram_run_fun = None
ans = 0
el = None
l_ans = None
n_list_v = list(oxi_info)
d = 0

def run_fun(name):
    global n_list_v
    global ans
    global el
    global l_ans
    n_list_v = []
    def find_ans():
        global ans
        global el
        global l_ans
        if el != None:
            el.destroy()
        x = None
        e = '('
        for i in range(len(v[1])):
            e += ev[i].get() + ','
        e = e[:-1:]
        e += ')'
        for i in oxi_info:
            if i[0] == name:
                x = i[1]
                break
        ans = eval('f_p.'+x[:x.index('(')]+e)[1]
        if l_ans != None:
            l_ans.destroy()
        l_ans = Label(fram_run_fun,text = f"Answer is {v[0][:v[0].index('=')+1]}{ans}")
        l_ans.grid(row = u+2,column = 0)
    def check_run_fun():
        a = True
        z = 0
        global el
        global l_ans
        for i in ev:
            if len(i.get()) == 0:
                #print(i.get)
                #el = Label(fram_run_fun,text = 'fill all values')
                messagebox.showinfo('entry sugession','fill all values')
                if l_ans != None:
                    l_ans.destroy()
                #el.grid(row = u+2,column = 0)
                #print(i.get()=='')
                a = False
                break
            #elif i.get().isdigit():
            #    #print('in')
            #    a = True
            find_ans()

    global fram_for
    global fram_run_fun
    fram_for.destroy()
    fram_run_fun = Frame(win)
    fram_run_fun.grid(row = 0,column = 0)
    for i in oxi_info:
        if i[0] == name:
            #print(i[1])
            v = oxi.run_and_return_function_formula_var(i[1])
            #print(v[1])
            break
    g = Label(fram_run_fun,text = 'Formula '+v[0])
    g.grid(row = 0,column= 0)
    ev = []
    s = 0
    u = 2
    eva = {}
    #print(v)
    for i in range(len(v[1])):
        ev.append(Entry(fram_run_fun))
    for i in ev:
        Label(fram_run_fun,text = f"enter the value of the veriable {v[1][s]}").grid(row = u-1,column = 0)
        i.grid(row = u,column = 0)
        u += 2
        s += 1
    k = Button(fram_run_fun,text = 'calculate',command = check_run_fun)
    k.grid(row = u+1,column= 0)
    #Label(fram_run_fun,text = f"Answer is {ans}").grid(row = len(v[1])+2,column = 0)
    h = Button(fram_run_fun,text = 'press to go back formula page',command = view_for)
    h.grid(row = u+4,column= 0)

def intro():
    global fram_intro
    global fram_for
    global fram_gra_con
    global fram_
    global fram_run_fun

    if fram_for != None:
        fram_for.destroy()
    if fram_gra_con != None:
        fram_gra_con.destroy()
    if fram_ != None:
        fram_.destroy()

    fram_intro = Frame(win)
    fram_intro.grid(row = 0,column = 0)
    Label(fram_intro,text = 'physics numerical solver',font = ('Times New Roman',15)).grid(row = 0,column = 0)
    Label(fram_intro,text = '',pady = 125).grid(row = 1,column = 0)
    but_for = Button(fram_intro,text = 'view formula to solve',command = view_for,font=('Times New Roman',17),bg = 'light blue')
    but_for.grid(row = 2,column = 0,padx = 40,pady = 10)
    but_gr_ap = Button(fram_intro,text = 'graph and real life applications',command = view_gra_con,font=('Times New Roman',17),bg = 'light blue')
    but_gr_ap.grid(row = 3,column = 0,padx = 40,pady = 10)
    but_hello = Button(fram_intro,text = '<not set up to now>',command = view_,font=('Times New Roman',17),bg = 'light blue')
    but_hello.grid(row = 4,column = 0,padx = 40,pady = 10)

def view_for():
    global n_list_v
    def abc():
        if list_box.get(ANCHOR) != '':
            run_fun(list_box.get(ANCHOR))
        else:
            Label(fram_for,text = 'select any forrmula to work',pady = 10).grid(row = 5,column = 0)

    def sca():
        global n_list_v
        n_list_v = []
        nm = 0
        #print('hello')
        for i in oxi_info:
            list_box.delete(0,END)
            if len(s_e.get()) != 0:
                if s_e.get() in i[0]:
                    n_list_v.append((i[0]))
                    #list_box.insert(END,i[0])
            else:
                n_list_v.append((i[0]))
        #print(n_list_v)
        #print(type(n_list_v))
        for i in n_list_v:
            list_box.insert(END,i)


    global v
    global fram_for
    global fram_run_fun
    v = ''
    fram_intro.destroy()
    if fram_run_fun != None:
        fram_run_fun.destroy()
    fram_for = Frame(win)
    fram_for.grid(row = 0,column = 0)
    Label(fram_for,text = 'in view formula',pady = 0).grid(row = 0,column = 0)
    Button(fram_for,text = 'press me to go back',command = intro).grid(row = 1,column = 0)
    Label(fram_for,text = 'search formula ').grid(row = 2,column = 0)
    s_e = Entry(fram_for)
    s_e.grid(row = 2,column = 1)
    Button(fram_for,text = 'search',command = sca).grid(row = 2,column = 3)
    list_box = Listbox(fram_for,width = 60)
    list_box.grid(row = 3,column = 0,columnspan = 4)
    #print(n_list_v)
    if len(n_list_v) == 0:
        n_list_v = oxi_info
    for i in n_list_v:
        list_box.insert(END,i[0])
    Button(fram_for,text = 'select',command = abc).grid(row = 4,column = 0)

def graph(f):
    global d
    def abcdefghij():
        print('hello')
        global d
        s = '('
        for i in range(len(vg[0])):
            s += str(el[i].get()) + ','
        s = s[:-1]
        s += ')'
        d = d.replace(d[d.index('('):],s)
        print(d)
        eval('g_p.'+d)
    fram_gra_con.destroy()
    fram_graph = Frame(win)
    fram_graph.grid(row = 0,column =0)
    vg = list()
    el = list()
    r = 0
    d = 0
    for i in oxi_gra_info:
        if i[0] == f:
            d = i[1]
            vg.append(i[1][i[1].index('(')+1:i[1].index(')')].split(','))
            vg.append(i)
    print(vg,end = '\n\n\n')
    for i in range(len(vg[0])):
        Label(fram_graph,text = f"enter the value of {vg[0][i]}").grid(row = r,column = 0)
        a = Entry(fram_graph)
        a.grid(row = r+1,column = 0)
        el.append(a)
        r += 2
    print('ebjs')
    Button(fram_graph,text = 'create graph',command = abcdefghij).grid(row = r,column = 0)

def view_gra_con():
    global fram_gra_con
    fram_intro.destroy()
    fram_gra_con = Frame(win)
    fram_gra_con.grid(row = 0,column = 0)
    Label(fram_gra_con,text = 'in view graph',pady = 0).grid(row = 0,column = 0)
    list_box_gra = Listbox(fram_gra_con)
    list_box_gra.grid(row = 1,column = 0)
    for i in oxi_gra_info:
        list_box_gra.insert(END,i[0])
    Button(fram_gra_con,text = 'view graph',command = lambda :graph(list_box_gra.get(ANCHOR))).grid(row = 2,column = 0)
    Button(fram_gra_con,text = 'press me to go back',command = intro).grid(row = 3,column = 0)

def view_():
    global fram_
    fram_intro.destroy()
    fram_ = Frame(win)
    fram_.grid(row = 0,column = 0)
    Label(fram_,text = 'in view -',pady = 0).grid(row = 0,column = 0)
    Button(fram_,text = 'press me to go back',command = intro).grid(row = 1,column = 0)


intro()
win.mainloop()
