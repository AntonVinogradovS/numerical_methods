import sys
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QPushButton,QTableWidget, QTableWidgetItem, QMainWindow, QLineEdit,QGraphicsScene
from PyQt5.QtGui import QPixmap
import numpy as np
import sympy
from sympy import *
from sympy.abc import x, y
from sympy.matrices import eye
from math import sin, cos, pi, e
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from _collections_abc import Iterable, Callable
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib
matplotlib.use('Qt5Agg')
sympy.init_printing()


#class MplCanvas(FigureCanvasQTAgg):
#
#    def __init__(w, parent=None, width=5, height=4, dpi=100):
#        fig = Figure(figsize=(width, height), dpi=dpi)
#        w.axes = fig.add_subplot(111)
#        super(MplCanvas, w).__init__(fig)
#
#class MplCanvas1(FigureCanvasQTAgg):
#
#    def __init__(w, parent=None, width=5, height=4, dpi=100):
#        fig = Figure(figsize=(width, height), dpi=dpi)
#        w.axes = fig.add_subplot(111)
#        super(MplCanvas, w).__init__(fig)


def calculator(n, h):
    res = []
    # Поиск численного решения
    A = sympy.Matrix([[-500.005, 499.995],[499.995, -500.005]])
    x = []
    v = []
    E = eye(2)
    x.append(0)
    tmp  = sympy.Matrix([7,13])
    v.append(tmp)
    for i in range(1, n):
        tempX = round(x[i - 1] + h, len(str(h))-2)
        x.append(tempX)
        vect = v[i-1]
        k1 = A*vect
        k2 = np.linalg.inv(np.asarray(E - A*h*0.5).astype('float64'))*(E + A*h*0.5)*A*vect
        t = vect + (k1 + k2)*h*0.5   
        v.append(t)
    # Поиск точного решения
    u = []
    A = sympy.Matrix(A)
    temp = A.eigenvects()
    alfa1 = 250000000000000/58925565098879
    alfa2 = -2500000000000000/176776695296637
    lambda1 = temp[0][0]
    lambda2 = temp[1][0]
    W1 = temp[0][2][0]
    W2 = temp[1][2][0]
    for i in range(n):
        t = alfa1*W1*e**(lambda1*x[i]) + alfa2*W2*e**(lambda2*x[i])
        u.append(t)
    # Расчет глобальной погрешности
    E_1 = []
    E_2 = []
    for i in range(n):
       E_1.append(u[i][0] - v[i][0])
       E_2.append(u[i][1] - v[i][1])
    # Создадим кортежи для удобной визуализации
    for i in range(n):
       currTuple = (i, x[i], u[i][0], u[i][1], v[i][0], v[i][1], E_1[i], E_2[i])
       res.append(currTuple)
    time = np.linspace(0, max(x), n)
    U1 = []
    U2 = []
    V1 = []
    V2 = []
    for i in range(n):
        U1.append(u[i][0])
        U2.append(u[i][1])
        V1.append(v[i][0])
        V2.append(v[i][1])  
    return res, U1, U2, V1, V2,x


def bildTable():
            # and one row
    
    # Set the table headers
    n = InputOne.text()
    n = int(n)
    h = InputTwo.text()
    h = float(h)
    table.setRowCount(n)
    res, u1, u2, v1, v2, x = calculator(n, h)
    #Set the tooltips to headings
    table.horizontalHeaderItem(0).setToolTip("Column 1 ")
    table.horizontalHeaderItem(1).setToolTip("Column 2 ")
    table.horizontalHeaderItem(2).setToolTip("Column 3 ")
    for i in range(n):
        table.setItem(i, 0, QTableWidgetItem(str(res[i][0])))
        table.setItem(i, 1, QTableWidgetItem(str(res[i][1])))    
        table.setItem(i, 2, QTableWidgetItem(str(res[i][2])))    
        table.setItem(i, 3, QTableWidgetItem(str(res[i][3])))    
        table.setItem(i, 4, QTableWidgetItem(str(res[i][4])))    
        table.setItem(i, 5, QTableWidgetItem(str(res[i][5])))    
        table.setItem(i, 6, QTableWidgetItem(str(res[i][6])))    
        table.setItem(i, 7, QTableWidgetItem(str(res[i][7])))
    
    time = np.linspace(0, max(x), n)
    plt.figure(figsize=(4, 4))
    plt.title("Точное решение")
    plt.plot(time, u1)
    plt.plot(time, u2)
    plt.xlabel('x')
    plt.ylabel('U1(x) и U2(x)')
    plt.savefig('saved_figure.png')
    plt.clf()

    plt.plot(time, v1)
    plt.plot(time, v2)
    plt.title("Численное решение")
    plt.xlabel('x')
    plt.ylabel('V1(x) и V2(x)')
    plt.savefig('saved_figure1.png')
    plt.clf()
    label = QLabel(w)
    pixmap = QPixmap("saved_figure.png")
    label.setPixmap(pixmap)
    label.move(10,180)
    
    label.show()

    label1 = QLabel(w)
    pixmap = QPixmap("saved_figure1.png")
    label1.setPixmap(pixmap)
    label1.move(500,180)
    label1.show()
    #fig, ax1 = plt.subplots()
    #data = np.array([0.7, 0.7, 0.7, 0.8, 0.9, 0.9, 1.5, 1.5, 1.5, 1.5])
    #bins = np.arange(0.6, 1.62, 0.02)
    #n1, bins1, patches1 = ax1.hist(data, bins, alpha=0.6, density=False, cumulative=False)
    # plot
    #w.plotWidget = FigureCanvas(fig)
    #w.setCentralWidget(w.plotWidget)
    #





#
    #sc2 = MplCanvas1(w, width=5, height=4, dpi=100)
    #sc2.axes.plot(x, V1)
    #sc2.axes.plot(x, V2)
    #sc2.show()
    #w.setCentralWidget(sc)
    #table.setItem(0, 0, QTableWidgetItem("Text in column 1"))
    #table.setItem(0, 1, QTableWidgetItem("Text in column 2"))
    #table.setItem(0, 2, QTableWidgetItem("Text in column 3"))

#def bildPlot():
    


    
 
if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = QWidget()
    w.resize(800,800)
    w.setWindowTitle('Жёсткая система дифференциальных уравнений')
    
    label = QLabel(w)
    label.setText("Кликай!")
    label.move(130,130)
    label.show()

# Ввод кол-ва шаговв
    labelCountStep = QLabel(w)
    labelCountStep.setText("Введите количество шагов: ")
    labelCountStep.move(10,10)
    labelCountStep.resize(150,10)
    labelCountStep.show()

    InputOne = w.textbox = QLineEdit(w)
    w.textbox.move(170, 12)
    w.textbox.resize(30,10)
# Ввод шага
    labelStep = QLabel(w)
    labelStep.setText("Введите размер шага: ")
    labelStep.move(10,35)
    labelStep.resize(140,10)
    labelStep.show()


    InputTwo = w.textbox = QLineEdit(w)
    w.textbox.move(170, 35)
    w.textbox.resize(30,10)
   
    #text, okPressed = QInputDialog.getText(w, "Get text","Your name:", QLineEdit.Normal, "")


    

    btn = QPushButton(w)
    btn.setText('Построить')
    btn.move(110,150)
    btn.show()
    btn.clicked.connect(bildTable)
    

    table = QTableWidget(w)
    table.move(1000,50)
    table.resize(850,500)
    table.setColumnCount(8)
    table.setHorizontalHeaderLabels(["n", "x", "U1", "U2", "V1", "V2", "E1", "E2"])
    
    

    

    #btnPlot = QPushButton(w)
    #btnPlot.setText('Построить график')
    #btnPlot.move(100,800)
    #btnPlot.show()
    #btnPlot.clicked.connect(bildPlot)

    

    

    
    
    
    w.show()
    sys.exit(app.exec_())