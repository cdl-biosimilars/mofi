# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ImportCsv.ui'
#
# Created by: PyQt5 UI code generator 5.9
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ImportCsvDialog(object):
    def setupUi(self, ImportCsvDialog):
        ImportCsvDialog.setObjectName("ImportCsvDialog")
        ImportCsvDialog.resize(384, 617)
        ImportCsvDialog.setModal(True)
        self.verticalLayout = QtWidgets.QVBoxLayout(ImportCsvDialog)
        self.verticalLayout.setSpacing(12)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox = QtWidgets.QGroupBox(ImportCsvDialog)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setMinimumSize(QtCore.QSize(0, 150))
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.teFileContents = QtWidgets.QTextEdit(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.teFileContents.sizePolicy().hasHeightForWidth())
        self.teFileContents.setSizePolicy(sizePolicy)
        self.teFileContents.setObjectName("teFileContents")
        self.verticalLayout_3.addWidget(self.teFileContents)
        self.verticalLayout.addWidget(self.groupBox)
        self.gridGroupBox = QtWidgets.QGroupBox(ImportCsvDialog)
        self.gridGroupBox.setObjectName("gridGroupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.gridGroupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.label_2 = QtWidgets.QLabel(self.gridGroupBox)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 3, 1, 1)
        self.sbSkipRows = QtWidgets.QSpinBox(self.gridGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sbSkipRows.sizePolicy().hasHeightForWidth())
        self.sbSkipRows.setSizePolicy(sizePolicy)
        self.sbSkipRows.setMinimumSize(QtCore.QSize(20, 0))
        self.sbSkipRows.setObjectName("sbSkipRows")
        self.gridLayout.addWidget(self.sbSkipRows, 2, 4, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(50, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 1, 2, 1, 1)
        self.leDecimal = QtWidgets.QLineEdit(self.gridGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leDecimal.sizePolicy().hasHeightForWidth())
        self.leDecimal.setSizePolicy(sizePolicy)
        self.leDecimal.setMinimumSize(QtCore.QSize(20, 0))
        self.leDecimal.setMaxLength(1)
        self.leDecimal.setObjectName("leDecimal")
        self.gridLayout.addWidget(self.leDecimal, 0, 4, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.gridGroupBox)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 0, 0, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.gridGroupBox)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 0, 3, 1, 1)
        self.leComment = QtWidgets.QLineEdit(self.gridGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leComment.sizePolicy().hasHeightForWidth())
        self.leComment.setSizePolicy(sizePolicy)
        self.leComment.setMinimumSize(QtCore.QSize(20, 0))
        self.leComment.setMaxLength(1)
        self.leComment.setObjectName("leComment")
        self.gridLayout.addWidget(self.leComment, 1, 1, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.gridGroupBox)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 1, 0, 1, 1)
        self.leSep = QtWidgets.QLineEdit(self.gridGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leSep.sizePolicy().hasHeightForWidth())
        self.leSep.setSizePolicy(sizePolicy)
        self.leSep.setMinimumSize(QtCore.QSize(20, 0))
        self.leSep.setMaxLength(1)
        self.leSep.setObjectName("leSep")
        self.gridLayout.addWidget(self.leSep, 0, 1, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.gridGroupBox)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 1, 3, 1, 1)
        self.leThousands = QtWidgets.QLineEdit(self.gridGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leThousands.sizePolicy().hasHeightForWidth())
        self.leThousands.setSizePolicy(sizePolicy)
        self.leThousands.setMinimumSize(QtCore.QSize(20, 0))
        self.leThousands.setMaxLength(1)
        self.leThousands.setObjectName("leThousands")
        self.gridLayout.addWidget(self.leThousands, 1, 4, 1, 1)
        self.cbHeader = QtWidgets.QCheckBox(self.gridGroupBox)
        self.cbHeader.setChecked(True)
        self.cbHeader.setObjectName("cbHeader")
        self.gridLayout.addWidget(self.cbHeader, 2, 0, 1, 2)
        spacerItem1 = QtWidgets.QSpacerItem(50, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem1, 1, 5, 1, 1)
        self.verticalLayout.addWidget(self.gridGroupBox)
        self.gbColumns = QtWidgets.QGroupBox(ImportCsvDialog)
        self.gbColumns.setMinimumSize(QtCore.QSize(0, 59))
        self.gbColumns.setObjectName("gbColumns")
        self.verticalLayout.addWidget(self.gbColumns)
        self.groupBox_2 = QtWidgets.QGroupBox(ImportCsvDialog)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.twPreview = QtWidgets.QTableWidget(self.groupBox_2)
        self.twPreview.setObjectName("twPreview")
        self.twPreview.setColumnCount(0)
        self.twPreview.setRowCount(0)
        self.verticalLayout_2.addWidget(self.twPreview)
        self.verticalLayout.addWidget(self.groupBox_2)
        self.buttonBox = QtWidgets.QDialogButtonBox(ImportCsvDialog)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Apply|QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(ImportCsvDialog)
        QtCore.QMetaObject.connectSlotsByName(ImportCsvDialog)

    def retranslateUi(self, ImportCsvDialog):
        _translate = QtCore.QCoreApplication.translate
        ImportCsvDialog.setWindowTitle(_translate("ImportCsvDialog", "Import CSV"))
        self.groupBox.setTitle(_translate("ImportCsvDialog", "Filename"))
        self.gridGroupBox.setTitle(_translate("ImportCsvDialog", "CSV options"))
        self.label_2.setText(_translate("ImportCsvDialog", "Skip rows:"))
        self.leDecimal.setText(_translate("ImportCsvDialog", "."))
        self.label_4.setText(_translate("ImportCsvDialog", "Delimiter:"))
        self.label_5.setText(_translate("ImportCsvDialog", "Decimal point:"))
        self.leComment.setText(_translate("ImportCsvDialog", "#"))
        self.label_3.setText(_translate("ImportCsvDialog", "Comment char:"))
        self.leSep.setText(_translate("ImportCsvDialog", ","))
        self.label_6.setText(_translate("ImportCsvDialog", "Thousands separator:"))
        self.cbHeader.setText(_translate("ImportCsvDialog", "Header"))
        self.gbColumns.setTitle(_translate("ImportCsvDialog", "Columns to use"))
        self.groupBox_2.setTitle(_translate("ImportCsvDialog", "Preview"))

