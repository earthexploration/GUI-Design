# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Bounds.ui'
#
# Created: Mon May 11 18:52:21 2015
#      by: PyQt4 UI code generator 4.9.6
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Bounds(QtGui.QWidget):
    def setupUi(self, Bounds):
        Bounds.setObjectName(_fromUtf8("Bounds"))
        Bounds.resize(838, 345)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Bounds.sizePolicy().hasHeightForWidth())
        Bounds.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setPointSize(10)
        Bounds.setFont(font)
        self.gridLayout_5 = QtGui.QGridLayout(Bounds)
        self.gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
        self.groupBox_Swift = QtGui.QGroupBox(Bounds)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_Swift.sizePolicy().hasHeightForWidth())
        self.groupBox_Swift.setSizePolicy(sizePolicy)
        self.groupBox_Swift.setObjectName(_fromUtf8("groupBox_Swift"))
        self.gridLayout = QtGui.QGridLayout(self.groupBox_Swift)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.lineEdit_C_Swift_U = QtGui.QLineEdit(self.groupBox_Swift)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_C_Swift_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_C_Swift_U.setSizePolicy(sizePolicy)
        self.lineEdit_C_Swift_U.setObjectName(_fromUtf8("lineEdit_C_Swift_U"))
        self.gridLayout.addWidget(self.lineEdit_C_Swift_U, 2, 2, 1, 1)
        self.lineEdit_A_Swift_U = QtGui.QLineEdit(self.groupBox_Swift)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_A_Swift_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_A_Swift_U.setSizePolicy(sizePolicy)
        self.lineEdit_A_Swift_U.setObjectName(_fromUtf8("lineEdit_A_Swift_U"))
        self.gridLayout.addWidget(self.lineEdit_A_Swift_U, 0, 2, 1, 1)
        self.label_B_Swift = QtGui.QLabel(self.groupBox_Swift)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_B_Swift.sizePolicy().hasHeightForWidth())
        self.label_B_Swift.setSizePolicy(sizePolicy)
        self.label_B_Swift.setObjectName(_fromUtf8("label_B_Swift"))
        self.gridLayout.addWidget(self.label_B_Swift, 1, 0, 1, 1)
        self.label_A_Swift = QtGui.QLabel(self.groupBox_Swift)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_A_Swift.sizePolicy().hasHeightForWidth())
        self.label_A_Swift.setSizePolicy(sizePolicy)
        self.label_A_Swift.setObjectName(_fromUtf8("label_A_Swift"))
        self.gridLayout.addWidget(self.label_A_Swift, 0, 0, 1, 1)
        self.lineEdit_B_Swift_U = QtGui.QLineEdit(self.groupBox_Swift)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_B_Swift_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_B_Swift_U.setSizePolicy(sizePolicy)
        self.lineEdit_B_Swift_U.setObjectName(_fromUtf8("lineEdit_B_Swift_U"))
        self.gridLayout.addWidget(self.lineEdit_B_Swift_U, 1, 2, 1, 1)
        self.lineEdit_A_Swift_L = QtGui.QLineEdit(self.groupBox_Swift)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_A_Swift_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_A_Swift_L.setSizePolicy(sizePolicy)
        self.lineEdit_A_Swift_L.setObjectName(_fromUtf8("lineEdit_A_Swift_L"))
        self.gridLayout.addWidget(self.lineEdit_A_Swift_L, 0, 1, 1, 1)
        self.label_C_Swift = QtGui.QLabel(self.groupBox_Swift)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_C_Swift.sizePolicy().hasHeightForWidth())
        self.label_C_Swift.setSizePolicy(sizePolicy)
        self.label_C_Swift.setObjectName(_fromUtf8("label_C_Swift"))
        self.gridLayout.addWidget(self.label_C_Swift, 2, 0, 1, 1)
        self.lineEdit_C_Swift_L = QtGui.QLineEdit(self.groupBox_Swift)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_C_Swift_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_C_Swift_L.setSizePolicy(sizePolicy)
        self.lineEdit_C_Swift_L.setObjectName(_fromUtf8("lineEdit_C_Swift_L"))
        self.gridLayout.addWidget(self.lineEdit_C_Swift_L, 2, 1, 1, 1)
        self.lineEdit_B_Swift_L = QtGui.QLineEdit(self.groupBox_Swift)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_B_Swift_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_B_Swift_L.setSizePolicy(sizePolicy)
        self.lineEdit_B_Swift_L.setObjectName(_fromUtf8("lineEdit_B_Swift_L"))
        self.gridLayout.addWidget(self.lineEdit_B_Swift_L, 1, 1, 1, 1)
        self.gridLayout_5.addWidget(self.groupBox_Swift, 0, 0, 1, 1)
        self.groupBox_Voce = QtGui.QGroupBox(Bounds)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_Voce.sizePolicy().hasHeightForWidth())
        self.groupBox_Voce.setSizePolicy(sizePolicy)
        self.groupBox_Voce.setObjectName(_fromUtf8("groupBox_Voce"))
        self.gridLayout_4 = QtGui.QGridLayout(self.groupBox_Voce)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.label_A_Voce = QtGui.QLabel(self.groupBox_Voce)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_A_Voce.sizePolicy().hasHeightForWidth())
        self.label_A_Voce.setSizePolicy(sizePolicy)
        self.label_A_Voce.setObjectName(_fromUtf8("label_A_Voce"))
        self.gridLayout_4.addWidget(self.label_A_Voce, 0, 0, 1, 1)
        self.lineEdit_A_Voce_L = QtGui.QLineEdit(self.groupBox_Voce)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_A_Voce_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_A_Voce_L.setSizePolicy(sizePolicy)
        self.lineEdit_A_Voce_L.setObjectName(_fromUtf8("lineEdit_A_Voce_L"))
        self.gridLayout_4.addWidget(self.lineEdit_A_Voce_L, 0, 1, 1, 1)
        self.lineEdit_A_Voce_U = QtGui.QLineEdit(self.groupBox_Voce)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_A_Voce_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_A_Voce_U.setSizePolicy(sizePolicy)
        self.lineEdit_A_Voce_U.setObjectName(_fromUtf8("lineEdit_A_Voce_U"))
        self.gridLayout_4.addWidget(self.lineEdit_A_Voce_U, 0, 2, 1, 1)
        self.label_B_Voce = QtGui.QLabel(self.groupBox_Voce)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_B_Voce.sizePolicy().hasHeightForWidth())
        self.label_B_Voce.setSizePolicy(sizePolicy)
        self.label_B_Voce.setObjectName(_fromUtf8("label_B_Voce"))
        self.gridLayout_4.addWidget(self.label_B_Voce, 1, 0, 1, 1)
        self.lineEdit_B_Voce_L = QtGui.QLineEdit(self.groupBox_Voce)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_B_Voce_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_B_Voce_L.setSizePolicy(sizePolicy)
        self.lineEdit_B_Voce_L.setObjectName(_fromUtf8("lineEdit_B_Voce_L"))
        self.gridLayout_4.addWidget(self.lineEdit_B_Voce_L, 1, 1, 1, 1)
        self.lineEdit_B_Voce_U = QtGui.QLineEdit(self.groupBox_Voce)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_B_Voce_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_B_Voce_U.setSizePolicy(sizePolicy)
        self.lineEdit_B_Voce_U.setObjectName(_fromUtf8("lineEdit_B_Voce_U"))
        self.gridLayout_4.addWidget(self.lineEdit_B_Voce_U, 1, 2, 1, 1)
        self.label_C_Voce = QtGui.QLabel(self.groupBox_Voce)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_C_Voce.sizePolicy().hasHeightForWidth())
        self.label_C_Voce.setSizePolicy(sizePolicy)
        self.label_C_Voce.setObjectName(_fromUtf8("label_C_Voce"))
        self.gridLayout_4.addWidget(self.label_C_Voce, 2, 0, 1, 1)
        self.lineEdit_C_Voce_L = QtGui.QLineEdit(self.groupBox_Voce)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_C_Voce_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_C_Voce_L.setSizePolicy(sizePolicy)
        self.lineEdit_C_Voce_L.setObjectName(_fromUtf8("lineEdit_C_Voce_L"))
        self.gridLayout_4.addWidget(self.lineEdit_C_Voce_L, 2, 1, 1, 1)
        self.lineEdit_C_Voce_U = QtGui.QLineEdit(self.groupBox_Voce)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_C_Voce_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_C_Voce_U.setSizePolicy(sizePolicy)
        self.lineEdit_C_Voce_U.setObjectName(_fromUtf8("lineEdit_C_Voce_U"))
        self.gridLayout_4.addWidget(self.lineEdit_C_Voce_U, 2, 2, 1, 1)
        self.gridLayout_5.addWidget(self.groupBox_Voce, 0, 1, 1, 1)
        self.groupBox_Gosh = QtGui.QGroupBox(Bounds)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_Gosh.sizePolicy().hasHeightForWidth())
        self.groupBox_Gosh.setSizePolicy(sizePolicy)
        self.groupBox_Gosh.setObjectName(_fromUtf8("groupBox_Gosh"))
        self.gridLayout_3 = QtGui.QGridLayout(self.groupBox_Gosh)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.label_A_Gosh = QtGui.QLabel(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_A_Gosh.sizePolicy().hasHeightForWidth())
        self.label_A_Gosh.setSizePolicy(sizePolicy)
        self.label_A_Gosh.setObjectName(_fromUtf8("label_A_Gosh"))
        self.gridLayout_3.addWidget(self.label_A_Gosh, 0, 0, 1, 1)
        self.lineEdit_A_Gosh_L = QtGui.QLineEdit(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_A_Gosh_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_A_Gosh_L.setSizePolicy(sizePolicy)
        self.lineEdit_A_Gosh_L.setObjectName(_fromUtf8("lineEdit_A_Gosh_L"))
        self.gridLayout_3.addWidget(self.lineEdit_A_Gosh_L, 0, 1, 1, 1)
        self.lineEdit_A_Gosh_U = QtGui.QLineEdit(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_A_Gosh_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_A_Gosh_U.setSizePolicy(sizePolicy)
        self.lineEdit_A_Gosh_U.setObjectName(_fromUtf8("lineEdit_A_Gosh_U"))
        self.gridLayout_3.addWidget(self.lineEdit_A_Gosh_U, 0, 2, 1, 1)
        self.label_B_Gosh = QtGui.QLabel(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_B_Gosh.sizePolicy().hasHeightForWidth())
        self.label_B_Gosh.setSizePolicy(sizePolicy)
        self.label_B_Gosh.setObjectName(_fromUtf8("label_B_Gosh"))
        self.gridLayout_3.addWidget(self.label_B_Gosh, 1, 0, 1, 1)
        self.lineEdit_B_Gosh_L = QtGui.QLineEdit(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_B_Gosh_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_B_Gosh_L.setSizePolicy(sizePolicy)
        self.lineEdit_B_Gosh_L.setObjectName(_fromUtf8("lineEdit_B_Gosh_L"))
        self.gridLayout_3.addWidget(self.lineEdit_B_Gosh_L, 1, 1, 1, 1)
        self.lineEdit_B_Gosh_U = QtGui.QLineEdit(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_B_Gosh_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_B_Gosh_U.setSizePolicy(sizePolicy)
        self.lineEdit_B_Gosh_U.setObjectName(_fromUtf8("lineEdit_B_Gosh_U"))
        self.gridLayout_3.addWidget(self.lineEdit_B_Gosh_U, 1, 2, 1, 1)
        self.label_C_Gosh = QtGui.QLabel(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_C_Gosh.sizePolicy().hasHeightForWidth())
        self.label_C_Gosh.setSizePolicy(sizePolicy)
        self.label_C_Gosh.setObjectName(_fromUtf8("label_C_Gosh"))
        self.gridLayout_3.addWidget(self.label_C_Gosh, 2, 0, 1, 1)
        self.lineEdit_C_Gosh_L = QtGui.QLineEdit(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_C_Gosh_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_C_Gosh_L.setSizePolicy(sizePolicy)
        self.lineEdit_C_Gosh_L.setObjectName(_fromUtf8("lineEdit_C_Gosh_L"))
        self.gridLayout_3.addWidget(self.lineEdit_C_Gosh_L, 2, 1, 1, 1)
        self.lineEdit_C_Gosh_U = QtGui.QLineEdit(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_C_Gosh_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_C_Gosh_U.setSizePolicy(sizePolicy)
        self.lineEdit_C_Gosh_U.setObjectName(_fromUtf8("lineEdit_C_Gosh_U"))
        self.gridLayout_3.addWidget(self.lineEdit_C_Gosh_U, 2, 2, 1, 1)
        self.label_D_Gosh = QtGui.QLabel(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_D_Gosh.sizePolicy().hasHeightForWidth())
        self.label_D_Gosh.setSizePolicy(sizePolicy)
        self.label_D_Gosh.setObjectName(_fromUtf8("label_D_Gosh"))
        self.gridLayout_3.addWidget(self.label_D_Gosh, 3, 0, 1, 1)
        self.lineEdit_D_Gosh_L = QtGui.QLineEdit(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_D_Gosh_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_D_Gosh_L.setSizePolicy(sizePolicy)
        self.lineEdit_D_Gosh_L.setObjectName(_fromUtf8("lineEdit_D_Gosh_L"))
        self.gridLayout_3.addWidget(self.lineEdit_D_Gosh_L, 3, 1, 1, 1)
        self.lineEdit_D_Gosh_U = QtGui.QLineEdit(self.groupBox_Gosh)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_D_Gosh_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_D_Gosh_U.setSizePolicy(sizePolicy)
        self.lineEdit_D_Gosh_U.setObjectName(_fromUtf8("lineEdit_D_Gosh_U"))
        self.gridLayout_3.addWidget(self.lineEdit_D_Gosh_U, 3, 2, 1, 1)
        self.gridLayout_5.addWidget(self.groupBox_Gosh, 1, 0, 1, 1)
        self.groupBox_HS = QtGui.QGroupBox(Bounds)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_HS.sizePolicy().hasHeightForWidth())
        self.groupBox_HS.setSizePolicy(sizePolicy)
        self.groupBox_HS.setObjectName(_fromUtf8("groupBox_HS"))
        self.gridLayout_2 = QtGui.QGridLayout(self.groupBox_HS)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.label_A_HS = QtGui.QLabel(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_A_HS.sizePolicy().hasHeightForWidth())
        self.label_A_HS.setSizePolicy(sizePolicy)
        self.label_A_HS.setObjectName(_fromUtf8("label_A_HS"))
        self.gridLayout_2.addWidget(self.label_A_HS, 0, 0, 1, 1)
        self.lineEdit_A_HS_L = QtGui.QLineEdit(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_A_HS_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_A_HS_L.setSizePolicy(sizePolicy)
        self.lineEdit_A_HS_L.setObjectName(_fromUtf8("lineEdit_A_HS_L"))
        self.gridLayout_2.addWidget(self.lineEdit_A_HS_L, 0, 1, 1, 1)
        self.lineEdit_A_HS_U = QtGui.QLineEdit(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_A_HS_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_A_HS_U.setSizePolicy(sizePolicy)
        self.lineEdit_A_HS_U.setObjectName(_fromUtf8("lineEdit_A_HS_U"))
        self.gridLayout_2.addWidget(self.lineEdit_A_HS_U, 0, 2, 1, 1)
        self.label_B_HS = QtGui.QLabel(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_B_HS.sizePolicy().hasHeightForWidth())
        self.label_B_HS.setSizePolicy(sizePolicy)
        self.label_B_HS.setObjectName(_fromUtf8("label_B_HS"))
        self.gridLayout_2.addWidget(self.label_B_HS, 1, 0, 1, 1)
        self.lineEdit_B_HS_L = QtGui.QLineEdit(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_B_HS_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_B_HS_L.setSizePolicy(sizePolicy)
        self.lineEdit_B_HS_L.setObjectName(_fromUtf8("lineEdit_B_HS_L"))
        self.gridLayout_2.addWidget(self.lineEdit_B_HS_L, 1, 1, 1, 1)
        self.lineEdit_B_HS_U = QtGui.QLineEdit(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_B_HS_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_B_HS_U.setSizePolicy(sizePolicy)
        self.lineEdit_B_HS_U.setObjectName(_fromUtf8("lineEdit_B_HS_U"))
        self.gridLayout_2.addWidget(self.lineEdit_B_HS_U, 1, 2, 1, 1)
        self.label_C_HS = QtGui.QLabel(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_C_HS.sizePolicy().hasHeightForWidth())
        self.label_C_HS.setSizePolicy(sizePolicy)
        self.label_C_HS.setObjectName(_fromUtf8("label_C_HS"))
        self.gridLayout_2.addWidget(self.label_C_HS, 2, 0, 1, 1)
        self.lineEdit_C_HS_L = QtGui.QLineEdit(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_C_HS_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_C_HS_L.setSizePolicy(sizePolicy)
        self.lineEdit_C_HS_L.setObjectName(_fromUtf8("lineEdit_C_HS_L"))
        self.gridLayout_2.addWidget(self.lineEdit_C_HS_L, 2, 1, 1, 1)
        self.lineEdit_C_HS_U = QtGui.QLineEdit(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_C_HS_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_C_HS_U.setSizePolicy(sizePolicy)
        self.lineEdit_C_HS_U.setObjectName(_fromUtf8("lineEdit_C_HS_U"))
        self.gridLayout_2.addWidget(self.lineEdit_C_HS_U, 2, 2, 1, 1)
        self.label_D_HS = QtGui.QLabel(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_D_HS.sizePolicy().hasHeightForWidth())
        self.label_D_HS.setSizePolicy(sizePolicy)
        self.label_D_HS.setObjectName(_fromUtf8("label_D_HS"))
        self.gridLayout_2.addWidget(self.label_D_HS, 3, 0, 1, 1)
        self.lineEdit_D_HS_L = QtGui.QLineEdit(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_D_HS_L.sizePolicy().hasHeightForWidth())
        self.lineEdit_D_HS_L.setSizePolicy(sizePolicy)
        self.lineEdit_D_HS_L.setObjectName(_fromUtf8("lineEdit_D_HS_L"))
        self.gridLayout_2.addWidget(self.lineEdit_D_HS_L, 3, 1, 1, 1)
        self.lineEdit_D_HS_U = QtGui.QLineEdit(self.groupBox_HS)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_D_HS_U.sizePolicy().hasHeightForWidth())
        self.lineEdit_D_HS_U.setSizePolicy(sizePolicy)
        self.lineEdit_D_HS_U.setObjectName(_fromUtf8("lineEdit_D_HS_U"))
        self.gridLayout_2.addWidget(self.lineEdit_D_HS_U, 3, 2, 1, 1)
        self.gridLayout_5.addWidget(self.groupBox_HS, 1, 1, 1, 1)

        self.retranslateUi(Bounds)
        QtCore.QMetaObject.connectSlotsByName(Bounds)
        Bounds.setTabOrder(self.lineEdit_A_Swift_L, self.lineEdit_A_Swift_U)
        Bounds.setTabOrder(self.lineEdit_A_Swift_U, self.lineEdit_B_Swift_L)
        Bounds.setTabOrder(self.lineEdit_B_Swift_L, self.lineEdit_B_Swift_U)
        Bounds.setTabOrder(self.lineEdit_B_Swift_U, self.lineEdit_C_Swift_L)
        Bounds.setTabOrder(self.lineEdit_C_Swift_L, self.lineEdit_C_Swift_U)
        Bounds.setTabOrder(self.lineEdit_C_Swift_U, self.lineEdit_A_Voce_L)
        Bounds.setTabOrder(self.lineEdit_A_Voce_L, self.lineEdit_A_Voce_U)
        Bounds.setTabOrder(self.lineEdit_A_Voce_U, self.lineEdit_B_Voce_L)
        Bounds.setTabOrder(self.lineEdit_B_Voce_L, self.lineEdit_B_Voce_U)
        Bounds.setTabOrder(self.lineEdit_B_Voce_U, self.lineEdit_C_Voce_L)
        Bounds.setTabOrder(self.lineEdit_C_Voce_L, self.lineEdit_C_Voce_U)
        Bounds.setTabOrder(self.lineEdit_C_Voce_U, self.lineEdit_A_Gosh_L)
        Bounds.setTabOrder(self.lineEdit_A_Gosh_L, self.lineEdit_A_Gosh_U)
        Bounds.setTabOrder(self.lineEdit_A_Gosh_U, self.lineEdit_B_Gosh_L)
        Bounds.setTabOrder(self.lineEdit_B_Gosh_L, self.lineEdit_B_Gosh_U)
        Bounds.setTabOrder(self.lineEdit_B_Gosh_U, self.lineEdit_C_Gosh_L)
        Bounds.setTabOrder(self.lineEdit_C_Gosh_L, self.lineEdit_C_Gosh_U)
        Bounds.setTabOrder(self.lineEdit_C_Gosh_U, self.lineEdit_D_Gosh_L)
        Bounds.setTabOrder(self.lineEdit_D_Gosh_L, self.lineEdit_D_Gosh_U)
        Bounds.setTabOrder(self.lineEdit_D_Gosh_U, self.lineEdit_A_HS_L)
        Bounds.setTabOrder(self.lineEdit_A_HS_L, self.lineEdit_A_HS_U)
        Bounds.setTabOrder(self.lineEdit_A_HS_U, self.lineEdit_B_HS_L)
        Bounds.setTabOrder(self.lineEdit_B_HS_L, self.lineEdit_B_HS_U)
        Bounds.setTabOrder(self.lineEdit_B_HS_U, self.lineEdit_C_HS_L)
        Bounds.setTabOrder(self.lineEdit_C_HS_L, self.lineEdit_C_HS_U)
        Bounds.setTabOrder(self.lineEdit_C_HS_U, self.lineEdit_D_HS_L)
        Bounds.setTabOrder(self.lineEdit_D_HS_L, self.lineEdit_D_HS_U)

    def retranslateUi(self, Bounds):
        Bounds.setWindowTitle(_translate("Bounds", "Bounds Setting", None))
        self.groupBox_Swift.setTitle(_translate("Bounds", "Swift : k(phi) = A*(B + phi)^C", None))
        self.lineEdit_C_Swift_U.setText(_translate("Bounds", "1", None))
        self.lineEdit_A_Swift_U.setText(_translate("Bounds", "1e+4", None))
        self.label_B_Swift.setText(_translate("Bounds", "B:", None))
        self.label_A_Swift.setText(_translate("Bounds", "A:", None))
        self.lineEdit_B_Swift_U.setText(_translate("Bounds", "1e-2", None))
        self.lineEdit_A_Swift_L.setText(_translate("Bounds", "0", None))
        self.label_C_Swift.setText(_translate("Bounds", "C:", None))
        self.lineEdit_C_Swift_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_B_Swift_L.setText(_translate("Bounds", "0", None))
        self.groupBox_Voce.setTitle(_translate("Bounds", "Voce : k(phi) = A + B*(1 - exp(-C*phi))", None))
        self.label_A_Voce.setText(_translate("Bounds", "A:", None))
        self.lineEdit_A_Voce_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_A_Voce_U.setText(_translate("Bounds", "1e+4", None))
        self.label_B_Voce.setText(_translate("Bounds", "B:", None))
        self.lineEdit_B_Voce_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_B_Voce_U.setText(_translate("Bounds", "1e+4", None))
        self.label_C_Voce.setText(_translate("Bounds", "C:", None))
        self.lineEdit_C_Voce_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_C_Voce_U.setText(_translate("Bounds", "1e+2", None))
        self.groupBox_Gosh.setTitle(_translate("Bounds", "Gosh : k(phi) = A*(B + phi)^C - D", None))
        self.label_A_Gosh.setText(_translate("Bounds", "A:", None))
        self.lineEdit_A_Gosh_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_A_Gosh_U.setText(_translate("Bounds", "1e+4", None))
        self.label_B_Gosh.setText(_translate("Bounds", "B:", None))
        self.lineEdit_B_Gosh_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_B_Gosh_U.setText(_translate("Bounds", "1e-2", None))
        self.label_C_Gosh.setText(_translate("Bounds", "C:", None))
        self.lineEdit_C_Gosh_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_C_Gosh_U.setText(_translate("Bounds", "1", None))
        self.label_D_Gosh.setText(_translate("Bounds", "D:", None))
        self.lineEdit_D_Gosh_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_D_Gosh_U.setText(_translate("Bounds", "1e+2", None))
        self.groupBox_HS.setTitle(_translate("Bounds", "Hockett_Sherby : k(phi) = A + B*(1 - exp(-C*phi^D))", None))
        self.label_A_HS.setText(_translate("Bounds", "A:", None))
        self.lineEdit_A_HS_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_A_HS_U.setText(_translate("Bounds", "1e+4", None))
        self.label_B_HS.setText(_translate("Bounds", "B:", None))
        self.lineEdit_B_HS_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_B_HS_U.setText(_translate("Bounds", "1e+4", None))
        self.label_C_HS.setText(_translate("Bounds", "C:", None))
        self.lineEdit_C_HS_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_C_HS_U.setText(_translate("Bounds", "1e+2", None))
        self.label_D_HS.setText(_translate("Bounds", "D:", None))
        self.lineEdit_D_HS_L.setText(_translate("Bounds", "0", None))
        self.lineEdit_D_HS_U.setText(_translate("Bounds", "1", None))
