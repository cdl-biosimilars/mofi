# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'importtabdata_ui.ui'
#
# Created by: PyQt5 UI code generator 5.9
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ImportTabData(object):
    def setupUi(self, ImportTabData):
        ImportTabData.setObjectName("ImportTabData")
        ImportTabData.resize(444, 646)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(ImportTabData.sizePolicy().hasHeightForWidth())
        ImportTabData.setSizePolicy(sizePolicy)
        ImportTabData.setModal(True)
        self.verticalLayout = QtWidgets.QVBoxLayout(ImportTabData)
        self.verticalLayout.setSpacing(12)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gbCsvFile = QtWidgets.QGroupBox(ImportTabData)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gbCsvFile.sizePolicy().hasHeightForWidth())
        self.gbCsvFile.setSizePolicy(sizePolicy)
        self.gbCsvFile.setMinimumSize(QtCore.QSize(0, 150))
        self.gbCsvFile.setObjectName("gbCsvFile")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.gbCsvFile)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.teFileContents = QtWidgets.QTextEdit(self.gbCsvFile)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.teFileContents.sizePolicy().hasHeightForWidth())
        self.teFileContents.setSizePolicy(sizePolicy)
        self.teFileContents.setMinimumSize(QtCore.QSize(0, 100))
        font = QtGui.QFont()
        font.setFamily("DejaVu Sans Mono")
        font.setPointSize(10)
        self.teFileContents.setFont(font)
        self.teFileContents.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.teFileContents.setReadOnly(True)
        self.teFileContents.setObjectName("teFileContents")
        self.verticalLayout_3.addWidget(self.teFileContents)
        self.verticalLayout.addWidget(self.gbCsvFile)
        self.gridGroupBox = QtWidgets.QGroupBox(ImportTabData)
        self.gridGroupBox.setObjectName("gridGroupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.gridGroupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.leSep = QtWidgets.QLineEdit(self.gridGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leSep.sizePolicy().hasHeightForWidth())
        self.leSep.setSizePolicy(sizePolicy)
        self.leSep.setMinimumSize(QtCore.QSize(20, 0))
        self.leSep.setObjectName("leSep")
        self.gridLayout.addWidget(self.leSep, 0, 1, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(50, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 0, 2, 3, 1)
        self.lbSkipRows = QtWidgets.QLabel(self.gridGroupBox)
        self.lbSkipRows.setObjectName("lbSkipRows")
        self.gridLayout.addWidget(self.lbSkipRows, 2, 3, 1, 1)
        self.sbSkipRows = QtWidgets.QSpinBox(self.gridGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sbSkipRows.sizePolicy().hasHeightForWidth())
        self.sbSkipRows.setSizePolicy(sizePolicy)
        self.sbSkipRows.setMinimumSize(QtCore.QSize(20, 0))
        self.sbSkipRows.setObjectName("sbSkipRows")
        self.gridLayout.addWidget(self.sbSkipRows, 2, 4, 1, 1)
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
        self.lbSep = QtWidgets.QLabel(self.gridGroupBox)
        self.lbSep.setObjectName("lbSep")
        self.gridLayout.addWidget(self.lbSep, 0, 0, 1, 1)
        self.lbDecimal = QtWidgets.QLabel(self.gridGroupBox)
        self.lbDecimal.setObjectName("lbDecimal")
        self.gridLayout.addWidget(self.lbDecimal, 0, 3, 1, 1)
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
        self.lbComment = QtWidgets.QLabel(self.gridGroupBox)
        self.lbComment.setObjectName("lbComment")
        self.gridLayout.addWidget(self.lbComment, 1, 0, 1, 1)
        self.lbThousands = QtWidgets.QLabel(self.gridGroupBox)
        self.lbThousands.setObjectName("lbThousands")
        self.gridLayout.addWidget(self.lbThousands, 1, 3, 1, 1)
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
        self.lbSheetName = QtWidgets.QLabel(self.gridGroupBox)
        self.lbSheetName.setObjectName("lbSheetName")
        self.gridLayout.addWidget(self.lbSheetName, 3, 3, 1, 1)
        self.cbSheetName = QtWidgets.QComboBox(self.gridGroupBox)
        self.cbSheetName.setObjectName("cbSheetName")
        self.gridLayout.addWidget(self.cbSheetName, 3, 4, 1, 1)
        self.cbHeader = QtWidgets.QCheckBox(self.gridGroupBox)
        self.cbHeader.setChecked(True)
        self.cbHeader.setObjectName("cbHeader")
        self.gridLayout.addWidget(self.cbHeader, 3, 0, 1, 1)
        self.lbQuote = QtWidgets.QLabel(self.gridGroupBox)
        self.lbQuote.setObjectName("lbQuote")
        self.gridLayout.addWidget(self.lbQuote, 2, 0, 1, 1)
        self.leQuote = QtWidgets.QLineEdit(self.gridGroupBox)
        self.leQuote.setMinimumSize(QtCore.QSize(20, 0))
        self.leQuote.setMaxLength(1)
        self.leQuote.setObjectName("leQuote")
        self.gridLayout.addWidget(self.leQuote, 2, 1, 1, 1)
        self.gridLayout.setColumnMinimumWidth(2, 50)
        self.verticalLayout.addWidget(self.gridGroupBox)
        self.groupBox_2 = QtWidgets.QGroupBox(ImportTabData)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.twPreview = QtWidgets.QTableWidget(self.groupBox_2)
        self.twPreview.setMinimumSize(QtCore.QSize(0, 100))
        self.twPreview.setObjectName("twPreview")
        self.twPreview.setColumnCount(0)
        self.twPreview.setRowCount(0)
        self.verticalLayout_2.addWidget(self.twPreview)
        self.verticalLayout.addWidget(self.groupBox_2)
        self.gbColumns = QtWidgets.QGroupBox(ImportTabData)
        self.gbColumns.setObjectName("gbColumns")
        self.verticalLayout.addWidget(self.gbColumns)
        self.buttonBox = QtWidgets.QDialogButtonBox(ImportTabData)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Apply|QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Help|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(ImportTabData)
        QtCore.QMetaObject.connectSlotsByName(ImportTabData)
        ImportTabData.setTabOrder(self.teFileContents, self.leSep)
        ImportTabData.setTabOrder(self.leSep, self.leComment)
        ImportTabData.setTabOrder(self.leComment, self.leQuote)
        ImportTabData.setTabOrder(self.leQuote, self.leDecimal)
        ImportTabData.setTabOrder(self.leDecimal, self.leThousands)
        ImportTabData.setTabOrder(self.leThousands, self.sbSkipRows)
        ImportTabData.setTabOrder(self.sbSkipRows, self.cbSheetName)
        ImportTabData.setTabOrder(self.cbSheetName, self.cbHeader)
        ImportTabData.setTabOrder(self.cbHeader, self.twPreview)

    def retranslateUi(self, ImportTabData):
        _translate = QtCore.QCoreApplication.translate
        ImportTabData.setWindowTitle(_translate("ImportTabData", "Import CSV"))
        self.gbCsvFile.setTitle(_translate("ImportTabData", "Filename"))
        self.gridGroupBox.setTitle(_translate("ImportTabData", "Options"))
        self.leSep.setText(_translate("ImportTabData", ","))
        self.lbSkipRows.setText(_translate("ImportTabData", "Skip rows:"))
        self.leDecimal.setText(_translate("ImportTabData", "."))
        self.lbSep.setText(_translate("ImportTabData", "Delimiter:"))
        self.lbDecimal.setText(_translate("ImportTabData", "Decimal point:"))
        self.leComment.setText(_translate("ImportTabData", "#"))
        self.lbComment.setText(_translate("ImportTabData", "Comment char:"))
        self.lbThousands.setText(_translate("ImportTabData", "Thousands separator:"))
        self.lbSheetName.setText(_translate("ImportTabData", "Sheet name:"))
        self.cbHeader.setText(_translate("ImportTabData", "Header"))
        self.lbQuote.setText(_translate("ImportTabData", "Quote char:"))
        self.leQuote.setText(_translate("ImportTabData", "\""))
        self.groupBox_2.setTitle(_translate("ImportTabData", "Preview"))
        self.gbColumns.setTitle(_translate("ImportTabData", "Columns to use"))

