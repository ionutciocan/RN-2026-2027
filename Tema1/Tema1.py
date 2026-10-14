import pathlib
import math


def load_system(path: pathlib.Path)->tuple[list[list[float]], list[float]]:
    matrix=[]
    B=[]
    with open(path,"r",encoding="utf-8") as f:
        for line in f:
            line=line.strip()
            if not line:
                continue
            left,right=line.split("=")
            l=[0.0,0.0,0.0]
            sign=1
            for term in left.split():
                if term == "+":
                    sign=1
                elif term == "-":
                    sign=-1
                else:
                    var=term[-1]
                    coef=term[:-1]
                    if coef in ("","+"):
                        val=1.0
                    elif coef=="-":
                        val=-1.0
                    else:
                        val=float(coef)
                    j="xyz".index(var)
                    l[j]+=sign*val
                    sign=1
            matrix.append(l)
            B.append(float(right))
        return matrix,B


def determinant(matrix: list[list[float]]) -> float:
    return matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1])-matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[1][2]*matrix[2][0])+matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[1][1]*matrix[2][0])

def trace(matrix: list[list[float]]) -> float:
    return matrix[0][0]+matrix[1][1]+matrix[2][2]

def norm(vector: list[float]) -> float:
    return math.sqrt(B[0]**2 + B[1]**2 + B[2]**2)

def transpose(matrix: list[list[float]]) -> list[list[float]]:
    return [
        [matrix[0][0], matrix[1][0], matrix[2][0]],
        [matrix[0][1], matrix[1][1], matrix[2][1]],
        [matrix[0][2], matrix[1][2], matrix[2][2]]
    ]
def multiply(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [matrix[0][0]*vector[0]+matrix[0][1]*vector[1]+matrix[0][2]*vector[2],matrix[1][0]*vector[0]+matrix[1][1]*vector[1]+matrix[1][2]*vector[2],matrix[2][0]*vector[0]+matrix[2][1]*vector[1]+matrix[2][2]*vector[2]]

def solve_cramer(matrix: list[list[float]], vector: list[float]) -> list[float]:
    dm=determinant(matrix)
    if dm==0:
        raise ValueError("Determinantul este zero; regula lui Cramer nu se aplică.")
    res=[]
    for c in range(0,3):
        cop=[line.copy() for line in matrix]
        for line in range(0,3):
            cop[line][c]=vector[line]
        res.append(determinant(cop)/dm)
    return res

def minor(matrix: list[list[float]], i: int, j: int) -> list[list[float]]:
    rezerv = [row.copy() for row in matrix]
    rezerv.pop(i)
    for row in rezerv:
        row.pop(j)
    return rezerv

def cofactor(matrix: list[list[float]]) -> list[list[float]]:
    res=list()
    for line in range(0,3):
        l=list()
        for col in range(0,3):
            sign=(-1)**(line+col)
            rezerv=minor(matrix,line,col)
            dt=rezerv[0][0]*rezerv[1][1]-rezerv[0][1]*rezerv[1][0]
            l.append(sign*dt)
        res.append(l)
    return res

def adjoint(matrix: list[list[float]]) -> list[list[float]]:
    trans=cofactor(matrix)
    return transpose(trans)

def solve(matrix: list[list[float]], vector: list[float]) -> list[float]:
    adj=adjoint(matrix)
    dt=determinant(matrix)
    for i in range(3):
        for j in range(3):
            adj[i][j]=adj[i][j]/dt
    return multiply(adj,vector)
    

A, B = load_system(pathlib.Path("system.txt"))
print(f"{A=} {B=}")
print(f"{determinant(A)=}")
print(f"{trace(A)=}")
print(f"{norm(B)=}")
print(f"{transpose(A)=}")
print(f"{multiply(A, B)=}")
print(f"{solve_cramer(A, B)=}")
print(f"{cofactor(A)=}")
print(f"{adjoint(A)=}")
print(f"{solve(A, B)=}")