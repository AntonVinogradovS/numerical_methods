import sys
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QPushButton,QTableWidget, QTableWidgetItem, QMainWindow, QLineEdit,QGraphicsScene, QTextEdit
from PyQt5.QtGui import QPixmap
import numpy as np
import sympy
from sympy import *
from sympy.abc import x, y
from sympy.matrices import eye
from math import  e
import math
import matplotlib.pyplot as plt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

import matplotlib
from IPython.display import display, Latex
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

def localError(x1, v1, h):
    A = sympy.Matrix([[-500.005, 499.995],[499.995, -500.005]])
    vv = v1
    E = eye(2)
    for i in range(1):
        vect = vv
        k1 = A*vect
        tt = (E + 0.5*h*A)#*A*vect
        t1 = E-h*0.5*A
        t1 = np.linalg.inv(np.asarray(t1).astype('float64'))
        k2 = t1*tt*A*vect
        #k2 = np.linalg.inv(np.asarray(E - A*h*0.5).astype('float64'))*(E + A*h*0.5)*A*vect
        t = vect + (k1 + k2)*h*0.5   
        #print(t)
    return t



def calculator(n, h, Correct):
    A = sympy.Matrix([[-500.005, 499.995],[499.995, -500.005]])
    x = []
    v = []
    hh = h
    E = eye(2)
    E
    x.append(0)
    tmp  = sympy.Matrix([7,13])
    v.append(tmp)
    eps = 0.001
    i = 1
    j = 0
    corr = []
    corr.append([0,0])
    flag = False
    res = []
    H = []
    H.append(0)
    multiplication = 0
    division = 0
    flag1 = False
    while i < n and math.isclose(v[i-1][0], Correct, rel_tol=1e-2) == False:
        #print(v[i-1])
        #print(h)

        vect = v[i-1]
        k1 = A*vect
        tt = (E + 0.5*h*A)#*A*vect
        t1 = E-h*0.5*A
        t1 = np.linalg.inv(np.asarray(t1).astype('float64'))
        k2 = t1*tt*A*vect
        #k2 = np.linalg.inv(np.asarray(E - A*h*0.5).astype('float64'))*(E + A*h*0.5)*A*vect
        t = vect + (k1 + k2)*h*0.5

        v_check = localError(x[-1], vect, h/2)
        v_check = localError(x[-1], v_check, h/2)

        S = (v_check - t)/3.


        S1 = S[0]
        S2 = S[1]
        #S = S[1]
        if vect[0] > 0.5 or flag1 == True:
            if eps/8 <= abs(S1) <= eps and eps/8 <= abs(S2) <= eps:
                #tempX = round(x[i - 1] + h, len(str(h))-2)
                tempX = x[i-1] + h
                x.append(tempX)
                #v.append(t)
                H.append(h)
                v.append(v_check)
                
                corr.append(t-v_check)
                
                
                i+=1
                continue
            elif abs(S1) < eps/8 and abs(S2) < eps/8:
                j +=1
                #tempX = round(x[i - 1] + h, len(str(h))-2)
                tempX = x[i-1] + h
                x.append(tempX)
                corr.append(t-v_check)
                #v.append(t)
                v.append(v_check)
                H.append(h)
                #print((t-v_check), h, i)
                print(v_check)
                h = h*2
                multiplication +=1
                i += 1
                continue
            else:
                division += 1
                h = h/2
                continue
        else:
            flag1 = True
            eps = 0.001
            h = h/2
            #tempX = round(x[i - 1] + h, len(str(h))-2)
            #x.append(tempX)
            #v.append(t)
            #H.append(h)
            #if flag == False:
            #    h = h/4
            #    division += 2
            #    flag = True
            #print(t , h)
            #v_check = localError(x[-1], t, h/2)
            #corr.append(t-v_check)
            #i+=1
        #поиск точного решения
    u = []
    A = sympy.Matrix(A)
    temp = A.eigenvects()
    alfa1 = 250000000000000/58925565098879
    alfa2 = -2500000000000000/176776695296637
    lambda1 = temp[0][0]
    lambda2 = temp[1][0]
    W1 = temp[0][2][0]
    W2 = temp[1][2][0]
    x1  = x
    n = i
    count = i
    for i in range(n):

        t = alfa1*W1*e**(lambda1*x1[i]) + alfa2*W2*e**(lambda2*x1[i])
        u.append(t)
    U1 = []
    U2 = []
    # Расчет глобальной погрешности
    E_1 = []
    E_2 = []
    for i in range(n):
       E_1.append(u[i][0] - v[i][0])
       E_2.append(u[i][1] - v[i][1])
    # Создадим кортежи для удобной визуализации
    for i in range(n):
       currTuple = (i, x[i], u[i][0], u[i][1], v[i][0], v[i][1], E_1[i], E_2[i], corr[i][0], corr[i][1])
       res.append(currTuple)
    max_e = 0
    min_e = 1
    max_E = 0
    min_E = -10
    for i in range(n):
        if abs(corr[i][0]) > max_e:
            max_e = abs(corr[i][0])
        if abs(corr[i][1]) > max_e:
            max_e = abs(corr[i][1])
        if abs(corr[i][0]) < min_e and corr[i][0] != 0:
            min_e = abs(corr[i][0])
        if abs(corr[i][1]) < min_e and corr[i][1] != 0:
            min_e = abs(corr[i][1])
        
        if abs(E_1[i]) > max_E:
            max_E = abs(E_1[i])
        if abs(E_2[i]) > max_E:
            max_E = abs(E_2[i])
        if abs(E_1[i]) < min_E:
            min_E = abs(E_1[i])
        if abs(E_2[i]) > min_E:
            min_E = abs(E_2[i])

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
    return res, U1, U2, V1, V2,x, count, H, multiplication, division, max_e, min_e, max_E, min_E









zzz = 0
def bildTable():
            # and one row
    
    # Set the table headers
    n = InputOne.text()
    n_MaxStep = n
    n = int(n)
    h = InputTwo.text()
    h_Firts = h
    h = float(h)
    
    global zzz
    
    zzz += 1
    Correct = InputFree.text()
    corStr = Correct
    Correct = float(Correct)
    
    res, u1, u2, v1, v2, x, count, H, multiplication, division, max_e, min_e, max_E, min_E = calculator(n, h, Correct)
    n = count


    
    
    labelMetod1.setText("Метод Рунге-Кутта неявный второго порядка " )
    labelMetod.setText("Входные данные :")
    labelX0.setText("Начальный шаг = " + h_Firts)
    labelHMax.setText("Максимальное количество шагов = " + n_MaxStep)
    labelEr.setText("Необходимая точность = " + corStr)
    labelRes.setText("Результат:")
    labelStep.setText("Количество шагов = " + str(count))
    labelStepMax.setText("Максильмальный шаг = " + str(max(H)))
    labelStepMin.setText("Минимальный шаг = " + str(min(H[1:])))
    ss = 'Последняя численная траектория: ' + "V1 = " + str(v1[-1]) + ", V2 = " + str(v2[-1]) + " при x = " + str(x[-1])
    labelEnd.setText(ss)
    labelGlobErr.setText("Максимальная / минимальная глобальная погрешность = " + str(max_E) + " / " + str(min_E))
    labelLocErr.setText("Максимальная / минимальная локальная погрешность = " + str(max_e) + " / " + str(min_e))
    labelMul.setText("Количество удвоений шага = " + str(multiplication))
    labelDiv.setText("Количество делений шага = " + str(division))




    table.setRowCount(n)
    #Set the tooltips to headings
    table.horizontalHeaderItem(0).setToolTip("Итерация")
    table.horizontalHeaderItem(1).setToolTip("Размер шага ")
    table.horizontalHeaderItem(2).setToolTip("Значение x")
    table.horizontalHeaderItem(3).setToolTip("Значение точной траектории (1)")
    table.horizontalHeaderItem(4).setToolTip("Значение точной траектории (2)")
    table.horizontalHeaderItem(5).setToolTip("Значение численной траектории (1)")
    table.horizontalHeaderItem(6).setToolTip("Значение численной траектории (2)")
    table.horizontalHeaderItem(7).setToolTip("Значение глобальной погрешности траектории (1)")
    table.horizontalHeaderItem(8).setToolTip("Значение глобальной погрешности траектории (2)")
    table.horizontalHeaderItem(9).setToolTip("Значение локальной погрешности траектории (1)")
    table.horizontalHeaderItem(10).setToolTip("Значение локальной погрешности траектории (2)")
    for i in range(n):
        table.setItem(i, 0, QTableWidgetItem(str(res[i][0])))

        table.setItem(i, 1, QTableWidgetItem(str(H[i])))

        table.setItem(i, 2, QTableWidgetItem(str(res[i][1])))    
        table.setItem(i, 3, QTableWidgetItem(str(res[i][2])))    
        table.setItem(i, 4, QTableWidgetItem(str(res[i][3])))    
        table.setItem(i, 5, QTableWidgetItem(str(res[i][4])))    
        table.setItem(i, 6, QTableWidgetItem(str(res[i][5])))    
        table.setItem(i, 7, QTableWidgetItem(str(res[i][6])))    
        table.setItem(i, 8, QTableWidgetItem(str(res[i][7])))
        table.setItem(i, 9, QTableWidgetItem(str(res[i][8])))    
        table.setItem(i, 10, QTableWidgetItem(str(res[i][9])))
    
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
    w.setWindowTitle('"Решение жёсткой системы ОДУ", Виноградов А.С., 382003_1')



    label = QLabel(w)
    label.setText("Кликай!")
    label.move(130,130)
    label.show()

    #w.init_handlers()

# Ввод кол-ва шаговв
    labelCountStep = QLabel(w)
    labelCountStep.setText("Введите количество шагов: ")
    labelCountStep.move(10,10)
    labelCountStep.resize(150,10)
    labelCountStep.show()

    InputOne = w.textbox = QLineEdit(w)
    w.textbox.move(170, 12)
    w.textbox.resize(40,20)
# Ввод шага
    labelStep = QLabel(w)
    labelStep.setText("Введите размер шага : ")
    labelStep.move(10,35)
    labelStep.resize(140,10)
    labelStep.show()
    labelStep.show()
    InputTwo = w.textbox = QLineEdit(w)
    w.textbox.move(170, 35)
    w.textbox.resize(40,20)


# Ввод шага
    labelCorrect = QLabel(w)
    labelCorrect.setText("Введите  точность: ")
    labelCorrect.move(10,60)
    labelCorrect.resize(150,10)
    labelCorrect.show()

    InputFree = w.textbox = QLineEdit(w)
    w.textbox.move(170, 58)
    w.textbox.resize(40,20)
   
    #text, okPressed = QInputDialog.getText(w, "Get text","Your name:", QLineEdit.Normal, "")

    labelPost = QLabel(w)
    labelPost.setText("Численно решить жёсткую систему ОДУ, используя неявный метод Рунге-Кутта второго порядка")
    labelPost.move(300,5)
    labelPost.resize(550,15)
    labelPost.show()

    label2 = QLabel(w)
    pixmap = QPixmap("C:/Users/rusan/AppData/Local/Microsoft/WindowsApps/PythonSoftwareFoundation.Python.3.10_qbz5n2kfra8p0/work/chm/render.png")
    label2.setPixmap(pixmap)
    label2.move(330,30)
    label2.show()

    label3 = QLabel(w)
    pixmap1 = QPixmap("C:/Users/rusan/AppData/Local/Microsoft/WindowsApps/PythonSoftwareFoundation.Python.3.10_qbz5n2kfra8p0/work/chm/render (1).png")
    label3.setPixmap(pixmap1)
    label3.move(340,110)
    label3.show()

    #справка
    

    

    btn = QPushButton(w)
    btn.setText('Построить')
    btn.move(110,150)
    btn.show()
    btn.clicked.connect(bildTable)
    
    

    table = QTableWidget(w)
    table.move(1000,50)
    table.resize(850,500)
    table.setColumnCount(11)
    table.setHorizontalHeaderLabels(["n", "h", "x", "U1", "U2", "V1", "V2", "E1", "E2", "e1", "e2"])
    
    labelMetod1 = QLabel(w)
    labelMetod1.move(1200,570)
    labelMetod1.resize(1000,10)
    labelMetod1.show()
    
    labelMetod = QLabel(w)
    labelMetod.move(1200,595)
    labelMetod.resize(1000,10)
    labelMetod.show()

    labelX0 = QLabel(w)
    labelX0.move(1200,625) 
    labelX0.resize(1000,10)
    labelX0.show()

    labelHMax = QLabel(w)
    labelHMax.move(1200,610)
    labelHMax.resize(1000,10)
    labelHMax.show()

    labelEr = QLabel(w)
    labelEr.move(1200,640)
    labelEr.resize(1000,10)
    labelEr.show()

    labelRes = QLabel(w)
    labelRes.move(1200,670)
    labelRes.resize(1000,10)
    labelRes.show()

    labelStep = QLabel(w)
    labelStep.move(1200,685)
    labelStep.resize(1000,10)
    labelStep.show()

    labelStepMax = QLabel(w)
    labelStepMax.move(1200,700)
    labelStepMax.resize(1000,10)
    labelStepMax.show()

    labelStepMin = QLabel(w)
    labelStepMin.move(1200,715)
    labelStepMin.resize(1000,10)
    labelStepMin.show()

    labelMul = QLabel(w)
    labelMul.move(1200, 730)
    labelMul.resize(1000,10)
    labelMul.show()

    labelDiv = QLabel(w)
    labelDiv.move(1200, 745)
    labelDiv.resize(1000,10)
    labelDiv.show()

    labelEnd = QLabel(w)
    labelEnd.move(1200,760)
    labelEnd.resize(1000,10)
    labelEnd.show()

    labelGlobErr = QLabel(w)
    labelGlobErr.move(1200,775)
    labelGlobErr.resize(1000,10)
    labelGlobErr.show()

    labelLocErr = QLabel(w)
    labelLocErr.move(1200,790)
    labelLocErr.resize(1000,10)
    labelLocErr.show()
    #btnPlot = QPushButton(w)
    #btnPlot.setText('Построить график')
    #btnPlot.move(100,800)
    #btnPlot.show()
    #btnPlot.clicked.connect(bildPlot)

    

    

    
    
    
    w.show()
    sys.exit(app.exec_())