# -*- coding: utf-8 -*-

import sys

from PyQt4 import QtCore,QtGui
from ui_mainwindow import Ui_MainWindow
from ui_bounds import Ui_Bounds
import pyqtgraph as pg

from reader import File
from flowcurves import FC
from yieldloci import YL
import numpy as np

import scipy.special._ufuncs_cxx

class BoundsBox(QtGui.QDialog):
    def __init__(self,parent=None):
        super(BoundsBox,self).__init__(parent)

        self.ui = Ui_Bounds(self)
        self.ui.setupUi(self)

class MainWindow(QtGui.QMainWindow):
    def __init__(self,parent=None):
        super(MainWindow,self).__init__(parent)

        self.ui = Ui_MainWindow(self)
        self.ui.setupUi(self)

        self.boundsbox = BoundsBox(self)
        self.boundsbox.hide()

        self.data = {} #data
        self.fn = None #filename
        self.xdata = None #temperature
        self.ydata = None #strainrate
        self.zdata = None #strain
        self.kdata = None #stress
        self.model = None #model
        self.fc = None #flowcurve
        self.params = None #parameters

        self.yl = None #yieldloci
        self.criterion = None #yieldcriterion
        
        self.ui.actionOpen.triggered.connect(self.openfile)
        self.ui.action.triggered.connect(self.about)
        self.connect(self.ui.actionExit,QtCore.SIGNAL('triggered()'),
                     QtCore.SLOT('close()'))

        self.ui.pushButton_FC.clicked.connect(self.on_pushButton_FC_clicked_)
        self.ui.pushButton_YL.clicked.connect(self.on_pushButton_YL_clicked_)
    def about(self):
        message = QtGui.QMessageBox.about(self,"About","This is a GUI of FlowCurve and Yield Loci in Sheet Metal Forming. --Hui Zhou")
        
    def openfile(self):
        filenames = QtGui.QFileDialog.getOpenFileNames(self, 'Open file',
                    '/home')
        for fn in filenames:
            try:
                data = File(fn)
                print "Import %s \n" %(fn)
                print "Check Temperature:"
                if data.temp_exp() is not None:
                    temp_n = len(data.temp_exp())
                    print temp_n
                else:
                    temp_n = 0
                    print None
                print "Check Strain Rate:"
                if data.rate_exp() is not None:
                    rate_n = len(data.rate_exp())
                    print rate_n
                else:
                    rate_n = 0
                    print None
                print "Check Strain:"
                if data.strainT_exp() is not None:
                    strain_n = len(data.strainT_exp())
                    print strain_n
                else:
                    strain_n = 0
                    print None
                print "Check Stress:"
                if data.stressT_exp() is not None:
                    stress_n = len(data.stressT_exp())
                    print stress_n
                else:
                    stress_n = 0
                    print None
                if strain_n !=0 and stress_n !=0 and strain_n == stress_n:
                    if temp_n == rate_n == 0:
                        self.data.setdefault(fn,data)
                        self.ui.comboBox_FN.addItem(fn)
                    if temp_n == rate_n == strain_n:
                        self.data.setdefault(fn,data)
                        self.ui.comboBox_FN.addItem(fn)
            except IOError:
                message = QtGui.QMessageBox.warning(self,"Warning","Read File:%s Failed!" %(fn))
                
    def on_comboBox_FN_activated(self):
        self.fn = self.ui.comboBox_FN.currentText()
        self.ui.comboBox_FN.setToolTip(self.fn)
        self.ui.comboBox_FN.setStatusTip(self.fn)
    def on_comboBox_FC_activated(self):
        self.model = str(self.ui.comboBox_FC.currentText())
        self.ui.pushButton_FC.setToolTip(self.model)
        self.ui.pushButton_FC.setStatusTip(self.model)
    def on_checkBox_Bounds_clicked(self):
        flag = self.ui.checkBox_Bounds.isChecked()
        if flag == True:
            self.boundsbox.show()
        if flag == False:
            self.boundsbox.hide()
    def on_checkBox_UE_clicked(self):
        self.on_comboBox_FN_activated()
        flag = self.ui.checkBox_UE.checkState()
        if flag == QtCore.Qt.Checked:
            try:
                self.xdata = self.data[self.fn].temp()
                self.ydata = self.data[self.fn].rate()
            except:
                self.xdata = None
                self.ydata = None
            self.zdata = self.data[self.fn].strainT()
            self.kdata = self.data[self.fn].stressT()
            return True
        if flag == QtCore.Qt.Unchecked:
            try:
                self.xdata = self.data[self.fn].temp_exp()
                self.ydata = self.data[self.fn].rate_exp()
            except:
                self.xdata = None
                self.ydata = None
            self.zdata = self.data[self.fn].strainT_exp()
            self.kdata = self.data[self.fn].stressT_exp()
            return False
    def on_pushButton_FC_clicked(self):
        pass
    def on_pushButton_FC_clicked_(self):
        
        self.ui.graphicsView_FC.clear()
        self.ui.graphicsView_FC.showGrid(True,True,alpha=0.7)
        self.on_checkBox_UE_clicked()
        self.on_comboBox_FC_activated()
        
        self.Line = pg.InfiniteLine(movable=True)
        self.Line.sigPositionChanged.connect(self.LineEdit)
        self.ui.graphicsView_FC.addItem(self.Line)
        self.ui.lineEdit_Xval.setText(str(self.Line.value()))
        self.ui.lineEdit_Xval.textChanged.connect(self.Xval)

        self.ui.graphicsView_FC.plot(self.zdata,self.kdata)
        try:
            bounds = {'Swift':[(float(self.boundsbox.ui.lineEdit_A_Swift_L.text()),float(self.boundsbox.ui.lineEdit_A_Swift_U.text())),
                               (float(self.boundsbox.ui.lineEdit_B_Swift_L.text()),float(self.boundsbox.ui.lineEdit_B_Swift_U.text())),
                               (float(self.boundsbox.ui.lineEdit_C_Swift_L.text()),float(self.boundsbox.ui.lineEdit_C_Swift_U.text()))],
                      'Voce':[(float(self.boundsbox.ui.lineEdit_A_Voce_L.text()),float(self.boundsbox.ui.lineEdit_A_Voce_U.text())),
                              (float(self.boundsbox.ui.lineEdit_B_Voce_L.text()),float(self.boundsbox.ui.lineEdit_B_Voce_U.text())),
                              (float(self.boundsbox.ui.lineEdit_C_Voce_L.text()),float(self.boundsbox.ui.lineEdit_C_Voce_U.text()))],
                      'Gosh':[(float(self.boundsbox.ui.lineEdit_A_Gosh_L.text()),float(self.boundsbox.ui.lineEdit_A_Gosh_U.text())),
                              (float(self.boundsbox.ui.lineEdit_B_Gosh_L.text()),float(self.boundsbox.ui.lineEdit_B_Gosh_U.text())),
                              (float(self.boundsbox.ui.lineEdit_C_Gosh_L.text()),float(self.boundsbox.ui.lineEdit_C_Gosh_U.text())),
                              (float(self.boundsbox.ui.lineEdit_D_Gosh_L.text()),float(self.boundsbox.ui.lineEdit_D_Gosh_U.text()))],
                      'Hockett_Sherby':[(float(self.boundsbox.ui.lineEdit_A_HS_L.text()),float(self.boundsbox.ui.lineEdit_A_HS_U.text())),
                                        (float(self.boundsbox.ui.lineEdit_B_HS_L.text()),float(self.boundsbox.ui.lineEdit_B_HS_U.text())),
                                        (float(self.boundsbox.ui.lineEdit_C_HS_L.text()),float(self.boundsbox.ui.lineEdit_C_HS_U.text())),
                                        (float(self.boundsbox.ui.lineEdit_D_HS_L.text()),float(self.boundsbox.ui.lineEdit_D_HS_U.text()))],
                      'SwiftVoce':[(float(self.boundsbox.ui.lineEdit_A_Swift_L.text()),float(self.boundsbox.ui.lineEdit_A_Swift_U.text())),
                                   (float(self.boundsbox.ui.lineEdit_B_Swift_L.text()),float(self.boundsbox.ui.lineEdit_B_Swift_U.text())),
                                   (float(self.boundsbox.ui.lineEdit_C_Swift_L.text()),float(self.boundsbox.ui.lineEdit_C_Swift_U.text())),
                                   (float(self.boundsbox.ui.lineEdit_A_Voce_L.text()),float(self.boundsbox.ui.lineEdit_A_Voce_U.text())),
                                   (float(self.boundsbox.ui.lineEdit_B_Voce_L.text()),float(self.boundsbox.ui.lineEdit_B_Voce_U.text())),
                                   (float(self.boundsbox.ui.lineEdit_C_Voce_L.text()),float(self.boundsbox.ui.lineEdit_C_Voce_U.text())),
                                   (0.0,1.0)],
                      }
                      
        except ValueError:
            message = QtGui.QMessageBox.warning(self,"Warning","Bounds Setting Type Error!")
        
        models = ["Swift","Voce","Gosh","Hockett_Sherby","SwiftVoce","Modified_Zener_Hollomon"]
        if self.model != u'Modified_Zener_Hollomon':
            fc = FC(self.zdata,self.kdata,ue=self.on_checkBox_UE_clicked)
            ans = fc.params(self.model,bounds[self.model])
            self.ui.label_FCErrorShow.setText(str(ans[1]))
            Y = fc.values(self.model,ans[0],self.zdata)
            y = fc.values(self.model,ans[0],float(self.Line.value()))
            self.ui.lineEdit_Yval.setText(str(y))
            self.ui.graphicsView_FC.plot(self.zdata,Y,pen='r')
            i = models.index(self.model)
            for j,p in enumerate(ans[0]):
                self.ui.treeWidget_FC.topLevelItem(i).child(j).setText(1,str(p))
        elif self.xdata and self.ydata is not None:
            fc = FC(self.zdata,self.kdata)
            ans = fc.Modified_Zener_HollomonFit(self.xdata,self.ydata,self.zdata,self.kdata)
            self.ui.label_FCErrorShow.setText(str(ans[1]))
            Y = fc.values('Modified_Zener_Hollomon',ans[0],self.zdata)
            y = fc.Modified_Zener_Hollomon(ans[0],sum(self.xdata)/len(self.xdata),sum(self.ydata)/len(self.ydata),float(self.Line.value()))
            self.ui.lineEdit_Yval.setText(str(y))
            self.ui.graphicsView_FC.plot(self.zdata,Y,pen='r')
            i = models.index(self.model)
            for j,p in enumerate(ans[0]):
                self.ui.treeWidget_FC.topLevelItem(i).child(j).setText(1,str(p))
        self.fc = fc
        self.params = ans[0]
        
    def LineEdit(self):
        self.ui.lineEdit_Xval.setText(str(self.Line.value()))
    def Xval(self):
        self.Line.setValue(float(self.ui.lineEdit_Xval.text()))
        if self.model != u'Modified_Zener_Hollomon':
            y = self.fc.values(self.model,self.params,float(self.ui.lineEdit_Xval.text()))
        else:
            y = self.fc.Modified_Zener_Hollomon(self.params,sum(self.xdata)/len(self.xdata),sum(self.ydata)/len(self.ydata),float(self.ui.lineEdit_Xval.text()))
        self.ui.lineEdit_Yval.setText(str(y))

    def on_comboBox_YL_activated(self):
        criterions = ["Vonmises","Hill48","Yld2000_2d","CPB"]
        self.criterion = str(self.ui.comboBox_YL.currentText())
        self.ui.tabWidget_Input.setCurrentIndex(criterions.index(self.criterion))
        self.ui.pushButton_YL.setToolTip(self.criterion)
        self.ui.pushButton_YL.setStatusTip(self.criterion)
    def on_pushButton_YL_clicked(self):
        pass
    def on_pushButton_YL_clicked_(self):
        
        self.ui.graphicsView_YL.clear()
        self.ui.graphicsView_YL.showGrid(True,True,alpha=0.7)
        self.ui.graphicsView_YL.setAspectLocked()
        self.on_comboBox_YL_activated()
        criterions = ["Vonmises","Hill48","Yld2000_2d","CPB"]
        try:
            Inputs = {"Vonmises": [float(self.ui.lineEdit_Vonmises_ay.text())],
                      "Hill48": [float(self.ui.lineEdit_Hill_a0.text()),
                                 float(self.ui.lineEdit_Hill_r0.text()),
                                 float(self.ui.lineEdit_Hill_r90.text())],
                      "Yld2000_2d": [float(self.ui.lineEdit_Yld_a0.text()),
                                     float(self.ui.lineEdit_Yld_a45.text()),
                                     float(self.ui.lineEdit_Yld_a90.text()),
                                     float(self.ui.lineEdit_Yld_ab.text()),
                                     float(self.ui.lineEdit_Yld_r0.text()),
                                     float(self.ui.lineEdit_Yld_r45.text()),
                                     float(self.ui.lineEdit_Yld_r90.text()),
                                     float(self.ui.lineEdit_Yld_rb.text()),
                                     float(self.ui.lineEdit_Yld_a.text())],
                      "CPB": [float(self.ui.lineEdit_CPB_a0_T.text()),
                              float(self.ui.lineEdit_CPB_a0_C.text()),
                              float(self.ui.lineEdit_CPB_a90_T.text()),
                              float(self.ui.lineEdit_CPB_a90_C.text()),
                              float(self.ui.lineEdit_CPB_aN_C.text()),
                              float(self.ui.lineEdit_CPB_ab_T.text()),
                              float(self.ui.lineEdit_CPB_ab_C.text()),
                              float(self.ui.lineEdit_CPB_a.text())],
                      }
        except ValueError:
            message = QtGui.QMessageBox.warning(self,"Warning","Inputs Setting Type Error!")

        yl = YL()
        if self.criterion == "Vonmises" and Inputs["Vonmises"] is not None:
            point = yl.Vonmises(Inputs["Vonmises"][0])
            self.ui.label_YLErrorShow.setText(str(0))
            coefficients = yl.params("Vonmises")
            self.ui.graphicsView_YL.plot(point['X'],point['Y'])
            
        if self.criterion == "Hill48" and Inputs["Hill48"] is not None:
            [a0,r0,r90] = Inputs["Hill48"]
            point = yl.Hill48(a0,None,r0,r90)
            self.ui.label_YLErrorShow.setText(str(0))
            coefficients = yl.params("Hill48")
            self.ui.graphicsView_YL.plot(point['X'],point['Y'])
            i = criterions.index(self.criterion)
            for j,p in enumerate(coefficients):
                self.ui.treeWidget_YL.topLevelItem(i).child(j).setText(1,str(p))
        if self.criterion == "Yld2000_2d" and Inputs["Yld2000_2d"] is not None:
            [a0,a45,a90,ab,r0,r45,r90,rb,a]=Inputs["Yld2000_2d"]
            point = yl.Yld2000_plot(a0,a45,a90,ab,r0,r45,r90,rb,a)
            coefficients = yl.params("Yld2000") #Yld2000
            err = yl.Yld2000_err(coefficients,a0,a45,a90,ab,r0,r45,r90,rb,a)
            error = np.sqrt(sum(np.square(err))/len(err))
            self.ui.label_YLErrorShow.setText(str(error))
            self.ui.graphicsView_YL.plot(point['X'],point['Y'])
            i = criterions.index(self.criterion)
            for j,p in enumerate(coefficients):
                self.ui.treeWidget_YL.topLevelItem(i).child(j).setText(1,str(p))
        if self.criterion == "CPB" and Inputs["CPB"] is not None:
            [a0_T,a0_C,a90_T,a90_C,aN_C,ab_T,ab_C,a] = Inputs["CPB"]
            point = yl.CPB_plot(a0_T,a0_C,a90_T,a90_C,aN_C,ab_T,ab_C,a)
            coefficients = yl.params("CPB")
            err = yl.CPB_err(coefficients,a0_T,a0_C,a90_T,a90_C,aN_C,ab_T,ab_C,a)
            error = np.sqrt(sum(np.square(err))/len(err))
            self.ui.label_YLErrorShow.setText(str(error))
            self.ui.graphicsView_YL.plot(point['X'],point['Y'])
            i = criterions.index(self.criterion)
            for j,p in enumerate(coefficients):
                self.ui.treeWidget_YL.topLevelItem(i).child(j).setText(1,str(p))
        
            
            


def main():
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec_()

if __name__ == '__main__':
    main()


        
