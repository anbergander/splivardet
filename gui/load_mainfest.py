import os

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QTableWidgetItem, QFileDialog, QMessageBox
from parameter import  Ui_Dialog


class Ui_MainWindow(object):
    def openWindow(self):
        for row in range(self.tableWidget.rowCount()):
            for column in range(self.tableWidget.columnCount()):
                item = self.tableWidget.item(row, column)
                if item is not None:
                    if column == 0:
                        if "barcode" in item.text():
                            continue
                        else:
                            QMessageBox.information(self.centralwidget, "Error", "Sample name must contain barcode + a number")
                            return
                    if column == 3:
                        file = item.text().split("\r")[0]
                        isdir = os.path.isfile(file)
                        if isdir:
                            continue
                        else:
                            QMessageBox.information(self.centralwidget, "Error", "Must be an existing path")
                            return
                else:
                    QMessageBox.information(self.centralwidget, "Error", "Table must be filled completely")
                    return
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_Dialog()
        self.ui.setupUi(self.window)
        self.ui.plainTextEdit_6.setPlainText(self.textEdit.toPlainText())
        self.window.show()

    def loadCsv(self):
        fileName, _ = QFileDialog.getOpenFileName(self.centralwidget, "Open TSV",
                                                  os.getenv('Home'), "TSV File (*.tsv) ;; TAB ( *.tab)")
        self.textEdit.setText(fileName)
        if fileName:
            self.loadCsvOnOpen(fileName)
    def loadCsvOnOpen(self, fileName):
        if fileName:
            df = pd.read_table(fileName, sep='\t', lineterminator="\n", header=None)
            print(df)
            header = df.iloc[0]
            self.tableWidget.setColumnCount(len(df.columns))
            self.tableWidget.setRowCount(len(df.index))
            for i in range(len(df.index)):
                for j in range(len(df.columns)):
                    self.tableWidget.setItem(i, j, QTableWidgetItem(str(df.iat[i, j])))
                    for j in range(len(df.columns)):
                        m = QTableWidgetItem(header[j])
                        self.tableWidget.setHorizontalHeaderItem(j, m)
            self.tableWidget.selectRow(0)
            self.tableWidget.resizeColumnsToContents()
            self.tableWidget.resizeRowsToContents()
            col_headers = ['Sample', 'Group', 'Batch', 'Path']
            self.tableWidget.setHorizontalHeaderLabels(col_headers)
            self.tableWidget.setColumnWidth(0, 110)
            self.tableWidget.setColumnWidth(1, 110)
            self.tableWidget.setColumnWidth(2, 110)
            self.tableWidget.setColumnWidth(3, 110)
            self.tableWidget.setDisabled(True)
    
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(540, 515)
        self.centralwidget = QtWidgets.QWidget(parent=MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.pushButton = QtWidgets.QPushButton(parent=self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(390, 450, 131, 41))
        self.pushButton.setObjectName("pushButton")
        font = QtGui.QFont()
        font.setPointSize(10)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("next-step-28.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.pushButton.setIcon(icon)
        self.pushButton.clicked.connect(self.openWindow)
        self.pushButton.clicked.connect(MainWindow.close)
        self.verticalLayoutWidget = QtWidgets.QWidget(parent=self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(30, 40, 471, 391))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.label_2 = QtWidgets.QLabel(parent=self.verticalLayoutWidget)
        font = QtGui.QFont()
        font.setPointSize(16)
        font.setFamily("Calibri")
        font.setBold(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_2.addWidget(self.label_2)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(parent=self.verticalLayoutWidget)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.textEdit = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget)
        self.textEdit.setObjectName("textEdit")
        self.textEdit.setMaximumHeight(30)
        self.textEdit.setDisabled(True)
        self.horizontalLayout.addWidget(self.textEdit)
        self.pushButton_3 = QtWidgets.QPushButton(parent=self.verticalLayoutWidget)
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.clicked.connect(self.loadCsv)
        self.horizontalLayout.addWidget(self.pushButton_3)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.tableWidget = QtWidgets.QTableWidget(parent=self.verticalLayoutWidget)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.verticalLayout.addWidget(self.tableWidget)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(parent=MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 540, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(parent=MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Load Manifest"))
        self.pushButton.setText(_translate("MainWindow", "Continue"))
        self.label_2.setText(_translate("MainWindow", "Load Manifest for Grouping"))
        self.label.setText(_translate("MainWindow", "Manifest"))
        self.textEdit.setPlaceholderText(_translate("MainWindow", "Path to Manifest"))
        self.pushButton_3.setText(_translate("MainWindow", "Browse"))

if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    Group = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(Group)
    Group.show()
    sys.exit(app.exec())