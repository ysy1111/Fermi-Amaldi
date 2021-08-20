import re

def parse_molecule (formula):
    pattern = r'[A-Z][a-z]?'
    pattern2 = r'\d+'
    map_dict = {'(': ')', '{': '}', '[': ']'}
    i = 0 
    temp_dict = {}
    while i < len(formula):
        ss = re.match(pattern, formula[i:])
        if ss is not None:
            if ss.group(0).isalpha():
                i += len(ss.group(0))
                ss2 = re.match(pattern2, formula[i:])
                if ss2 is not None:
                    if ss2.group(0).isdigit():
                        i += len(ss2.group(0))
                        temp_dict[ss.group(0)] = int(ss2.group(0)) + temp_dict.get(ss.group(0), 0)
                    else:
                        temp_dict[ss.group(0)] = 1 + temp_dict.get(ss.group(0), 0)
                else:       
                    temp_dict[ss.group(0)] = 1 + temp_dict.get(ss.group(0), 0)
        elif formula[i] in map_dict:
            index = formula[i:].index(map_dict[formula[i]])
            temp_dict2 = parse_molecule(formula[i+1:i+index])
            i += index + 1 
            ss2 = re.match(pattern2, formula[i:])
            if ss2 is not None:
                if ss2.group(0).isdigit():
                    i += len(ss2.group(0))
                    for item in temp_dict2:
                        temp_dict2[item] *= int(ss2.group(0))
            for item in temp_dict2:
                temp_dict[item] = temp_dict2[item] + temp_dict.get(item, 0)
    return temp_dict

def _element_list_(sys):
    if sys[-1]=='+':
        sys = sys[:-1]
    elif sys[-1]=='-':
        sys = sys[:-1]
    return parse_molecule(sys)


def _electron_number_(sys):
    Dict_for_nelec = {'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18}
    electron_number = 0
    if sys[-1]=='+':
        electron_number -=1
        sys = sys[:-1]
    elif sys[-1]=='-':
        electron_number +=1
        sys = sys[:-1]
    element = parse_molecule(sys)
    for i in element:
        electron_number += element[i]*Dict_for_nelec[i]
    return electron_number

def _atom_number_(sys):
    Dict_for_nelec = {'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18}
    atom_number = 0
    if sys[-1]=='-' or sys[-1]=='+':
        sys = sys[:-1]
    
    element = parse_molecule(sys)
    for i in element:
        atom_number += element[i]
    return atom_number