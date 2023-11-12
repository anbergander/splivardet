import csv
import os
import sys

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QFileDialog, QMessageBox
from de.de_final import main2


class Ui_Metadata(object):
    def browseDirectory(self):
        fname = QFileDialog.getExistingDirectory(self.centralwidget, "Select Directory", os.getenv('Home'))
        self.textEdit_2.setPlainText(fname)
    def browseCsv(self):
        fname = QFileDialog.getOpenFileName(self.centralwidget, "Open file", os.getenv('Home'), 'TSV (*.tsv) ;; TXT (*.txt) ;; CSV (*.csv)')
        self.textEdit.setPlainText(fname[0])

    def savetable(self):
        for row in range(self.tableWidget.rowCount()):
            for column in range(self.tableWidget.columnCount()):
                item = self.tableWidget.item(row, column)
                if item is not None:
                    continue
                else:
                    QMessageBox.information(self.centralwidget, "Error", "Table must be filled completely")
                    return
        path = self.textEdit_2.toPlainText() + '/metadata.tsv'
        print(path)
        if path != '':
            with open(path, 'w') as csv_file:
                writer = csv.writer(csv_file, delimiter="\t", lineterminator='\n')
                factors = self.textEdit_3.toPlainText()
                list = factors.split(",")
                list = ['Sample'] + list
                writer.writerow(list)
                for row in range(self.tableWidget.rowCount()):
                    row_data = []
                    for column in range(self.tableWidget.columnCount()):
                        item = self.tableWidget.item(row, column)
                        if item is not None:
                            row_data.append(item.text())
                        else:
                            row_data.append('')
                    writer.writerow(row_data)
    def checkifFilledBeforeContinue(self):
        if self.textEdit_2.toPlainText() == '' or self.textEdit.toPlainText() == '' or self.textEdit_3.toPlainText() == '':
            self.pushButton.setDisabled(True)
        else:
            self.pushButton.setDisabled(False)

    def checkIfFilled(self):
        if self.textEdit_3.toPlainText() == '':
            self.pushButton_4.setDisabled(True)
        else:
            self.pushButton_4.setDisabled(False)
    def performDe(self):
        self.pushButton.setDisabled(True)
        self.pushButton_2.setDisabled(True)
        self.pushButton_3.setDisabled(True)
        self.savetable()
        count = self.textEdit.toPlainText()
        metadata = self.textEdit_2.toPlainText() + '/metadata.tsv'
        main2(count, metadata, self.textEdit_2.toPlainText())
        QMessageBox.information(self.centralwidget, "Done", "DEG was performed")
    def createMetadataTable(self):
        self.tableWidget.show()
        factors = self.textEdit_3.toPlainText()
        self.textEdit_3.setDisabled(True)
        self.textEdit_3.hide()
        self.label_3.hide()
        self.spinBox.hide()
        self.label_4.hide()
        list = factors.split(",")
        list = ['Sample'] + list
        rows = self.spinBox.value()
        self.spinBox.setDisabled(True)
        self.pushButton_4.hide()
        self.tableWidget.setColumnCount(len(list))
        self.tableWidget.setRowCount(rows)
        self.tableWidget.setHorizontalHeaderLabels(list)

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.setFixedSize(816, 625)
        self.centralwidget = QtWidgets.QWidget(parent=MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.pushButton = QtWidgets.QPushButton(parent=self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(644, 550, 161, 31))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.setDisabled(True)
        self.pushButton.clicked.connect(self.performDe)
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(parent=self.centralwidget)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(40, 50, 731, 461))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_2 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_2)
        self.label_2.setObjectName("label_2")
        self.label_2.setMinimumWidth(100)
        self.horizontalLayout.addWidget(self.label_2)
        self.textEdit = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_2)
        self.textEdit.setObjectName("textBrowser")
        self.textEdit.setDisabled(True)
        self.textEdit.setMaximumHeight(50)
        self.textEdit.textChanged.connect(self.checkifFilledBeforeContinue)
        self.horizontalLayout.addWidget(self.textEdit)
        self.pushButton_2 = QtWidgets.QPushButton(parent=self.verticalLayoutWidget_2)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.clicked.connect(self.browseCsv)
        self.horizontalLayout.addWidget(self.pushButton_2)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label = QtWidgets.QLabel(parent=self.verticalLayoutWidget_2)
        self.label.setObjectName("label")
        self.label.setMinimumWidth(100)
        self.horizontalLayout_2.addWidget(self.label)
        self.textEdit_2 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_2)
        self.textEdit_2.setObjectName("textBrowser_2")
        self.textEdit_2.setDisabled(True)
        self.textEdit_2.textChanged.connect(self.checkifFilledBeforeContinue)
        self.textEdit_2.setMaximumHeight(50)
        self.horizontalLayout_2.addWidget(self.textEdit_2)
        self.pushButton_3 = QtWidgets.QPushButton(parent=self.verticalLayoutWidget_2)
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.clicked.connect(self.browseDirectory)
        self.horizontalLayout_2.addWidget(self.pushButton_3)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_3 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_2)
        self.label_3.setObjectName("label_3")
        self.label_3.setMinimumWidth(100)
        self.horizontalLayout_3.addWidget(self.label_3)
        self.textEdit_3 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_2)
        self.textEdit_3.setObjectName("lineEdit")
        self.textEdit_3.textChanged.connect(self.checkIfFilled)
        self.textEdit_3.textChanged.connect(self.checkifFilledBeforeContinue)
        self.horizontalLayout_3.addWidget(self.textEdit_3)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.label_4 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_2)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_4.addWidget(self.label_4)
        self.spinBox = QtWidgets.QSpinBox(parent=self.verticalLayoutWidget_2)
        self.spinBox.setObjectName("spinBox")
        self.spinBox.setMinimum(1)
        self.horizontalLayout_4.addWidget(self.spinBox)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        self.pushButton_4 = QtWidgets.QPushButton(parent=self.verticalLayoutWidget_2)
        self.pushButton_4.setObjectName("pushButton_4")
        self.pushButton_4.clicked.connect(self.createMetadataTable)
        self.pushButton_4.setDisabled(True)
        self.verticalLayout.addWidget(self.pushButton_4)
        self.tableWidget = QtWidgets.QTableWidget(parent=self.verticalLayoutWidget_2)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.tableWidget.hide()
        self.verticalLayout.addWidget(self.tableWidget)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.label_5 = QtWidgets.QLabel(parent=self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(250, 10, 300, 20))
        self.label_5.setObjectName("label_5")
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        self.label_5.setFont(font)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(parent=MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 816, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(parent=MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.pushButton.setText(_translate("MainWindow", "Start"))
        self.label_2.setText(_translate("MainWindow", "Count Table"))
        self.textEdit.setPlaceholderText(_translate("MainWindow", "Count Data"))
        self.pushButton_2.setText(_translate("MainWindow", "Browse"))
        self.label.setText(_translate("MainWindow", "Output"))
        self.textEdit_2.setPlaceholderText(_translate("MainWindow", "Output Path"))
        self.pushButton_3.setText(_translate("MainWindow", "Browse"))
        self.label_3.setText(_translate("MainWindow", "Factors"))
        self.textEdit_3.setPlaceholderText(_translate("MainWindow", "Factors for multivariate DEG"))
        self.label_4.setText(_translate("MainWindow", "Number of Samples"))
        self.pushButton_4.setText(_translate("MainWindow", "Generate Metadata Table"))
        self.label_5.setText(_translate("MainWindow", "Differential Gene Expression Analysis"))

if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    Window = QtWidgets.QMainWindow()
    ui = Ui_Metadata()
    ui.setupUi(Window)
    Window.show()
    sys.exit(app.exec())