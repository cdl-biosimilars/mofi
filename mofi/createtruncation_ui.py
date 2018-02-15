# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'createtruncation_ui.ui'
#
# Created by: PyQt5 UI code generator 5.9
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_CreateTruncation(object):
    def setupUi(self, CreateTruncation):
        CreateTruncation.setObjectName("CreateTruncation")
        CreateTruncation.resize(445, 333)
        self.verticalLayout = QtWidgets.QVBoxLayout(CreateTruncation)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(CreateTruncation)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.btFromParameters = QtWidgets.QPushButton(CreateTruncation)
        self.btFromParameters.setObjectName("btFromParameters")
        self.horizontalLayout.addWidget(self.btFromParameters)
        self.lbChain = QtWidgets.QLabel(CreateTruncation)
        self.lbChain.setObjectName("lbChain")
        self.horizontalLayout.addWidget(self.lbChain)
        self.cbChain = QtWidgets.QComboBox(CreateTruncation)
        self.cbChain.setObjectName("cbChain")
        self.horizontalLayout.addWidget(self.cbChain)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.teSequence = QtWidgets.QTextEdit(CreateTruncation)
        self.teSequence.setAcceptRichText(False)
        self.teSequence.setObjectName("teSequence")
        self.verticalLayout.addWidget(self.teSequence)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setVerticalSpacing(0)
        self.gridLayout.setObjectName("gridLayout")
        self.rbResidueNth = QtWidgets.QRadioButton(CreateTruncation)
        self.rbResidueNth.setChecked(True)
        self.rbResidueNth.setObjectName("rbResidueNth")
        self.bgResidue = QtWidgets.QButtonGroup(CreateTruncation)
        self.bgResidue.setObjectName("bgResidue")
        self.bgResidue.addButton(self.rbResidueNth)
        self.gridLayout.addWidget(self.rbResidueNth, 0, 3, 1, 1)
        self.rbNterminus = QtWidgets.QRadioButton(CreateTruncation)
        self.rbNterminus.setChecked(True)
        self.rbNterminus.setObjectName("rbNterminus")
        self.bgTerminus = QtWidgets.QButtonGroup(CreateTruncation)
        self.bgTerminus.setObjectName("bgTerminus")
        self.bgTerminus.addButton(self.rbNterminus)
        self.gridLayout.addWidget(self.rbNterminus, 0, 1, 1, 1)
        self.rbCterminus = QtWidgets.QRadioButton(CreateTruncation)
        self.rbCterminus.setObjectName("rbCterminus")
        self.bgTerminus.addButton(self.rbCterminus)
        self.gridLayout.addWidget(self.rbCterminus, 1, 1, 1, 1)
        self.rbResidueList = QtWidgets.QRadioButton(CreateTruncation)
        self.rbResidueList.setObjectName("rbResidueList")
        self.bgResidue.addButton(self.rbResidueList)
        self.gridLayout.addWidget(self.rbResidueList, 1, 3, 1, 1)
        self.leResidueList = QtWidgets.QLineEdit(CreateTruncation)
        self.leResidueList.setEnabled(False)
        self.leResidueList.setObjectName("leResidueList")
        self.gridLayout.addWidget(self.leResidueList, 1, 4, 1, 1)
        self.label_2 = QtWidgets.QLabel(CreateTruncation)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 0, 0, 2, 1)
        self.label_3 = QtWidgets.QLabel(CreateTruncation)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 0, 2, 2, 1)
        self.sbResidueNth = QtWidgets.QSpinBox(CreateTruncation)
        self.sbResidueNth.setMinimum(1)
        self.sbResidueNth.setMaximum(999)
        self.sbResidueNth.setProperty("value", 1)
        self.sbResidueNth.setObjectName("sbResidueNth")
        self.gridLayout.addWidget(self.sbResidueNth, 0, 4, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        spacerItem1 = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.verticalLayout.addItem(spacerItem1)
        self.buttonBox = QtWidgets.QDialogButtonBox(CreateTruncation)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(CreateTruncation)
        self.buttonBox.accepted.connect(CreateTruncation.accept)
        self.buttonBox.rejected.connect(CreateTruncation.reject)
        QtCore.QMetaObject.connectSlotsByName(CreateTruncation)

    def retranslateUi(self, CreateTruncation):
        _translate = QtCore.QCoreApplication.translate
        CreateTruncation.setWindowTitle(_translate("CreateTruncation", "Truncation"))
        self.label.setText(_translate("CreateTruncation", "Sequence:"))
        self.btFromParameters.setText(_translate("CreateTruncation", "From parameters"))
        self.lbChain.setText(_translate("CreateTruncation", "Chain:"))
        self.rbResidueNth.setText(_translate("CreateTruncation", "every n-th residue:"))
        self.rbNterminus.setText(_translate("CreateTruncation", "N"))
        self.rbCterminus.setText(_translate("CreateTruncation", "C"))
        self.rbResidueList.setText(_translate("CreateTruncation", "these residues:"))
        self.label_2.setText(_translate("CreateTruncation", "Truncate from"))
        self.label_3.setText(_translate("CreateTruncation", "terminus after"))

