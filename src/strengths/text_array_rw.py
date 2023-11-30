def load_1D_array_txt(path, elements_type) :
    f = open(path, "r")
    s = f.read()
    s = s.replace(",", " ")
    t = s.split()
    l = []
    for i in range(len(t)) : 
        l.append(elements_type(t[i]))
    return l

def save_1D_array_txt(l, path) :
    f = open(path, "w")
    s = ""
    for i in range(len(l)) : 
        s += str(l[i]) + " "
    f.write(s)