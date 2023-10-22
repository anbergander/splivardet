import os

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QFileDialog

from de import go_final


class Ui_GO(object):
    def checkIfFilled(self):
        if (self.textEdit.toPlainText() == '' or self.textEdit_2.toPlainText() == ''):
            self.pushButton_2.setDisabled(True)
        else:
            self.pushButton_2.setDisabled(False)
    def selectDeTable(self):
        fileName, _ = QFileDialog.getOpenFileName(self.verticalLayoutWidget, "Open CSV",
                                                  os.getenv('Home'), "CSV Files (*.csv);; TSV Files (*.tsv);; "
                                                                     "TXT Files (*.txt) ;; TAB Files (*.tab)")
        self.textEdit.setPlainText(fileName)
    def browseDirectoryForOutput(self):
        fname = QFileDialog.getExistingDirectory(self.verticalLayoutWidget, "Select Directory", os.getenv('Home'))
        self.textEdit_2.setPlainText(fname)
        print(fname)
    def performGO(self):
        go_final.main(self.textEdit_2.toPlainText(), self.textEdit.toPlainText())
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(470, 225)
        self.label = QtWidgets.QLabel(parent=Dialog)
        self.label.setGeometry(QtCore.QRect(200, 10, 100, 20))
        self.label.setObjectName("label")
        font = QtGui.QFont()
        font.setPointSize(16)
        font.setFamily("Calibri")
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.verticalLayoutWidget = QtWidgets.QWidget(parent=Dialog)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(20, 50, 430, 150))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_2 = QtWidgets.QLabel(parent=self.verticalLayoutWidget)
        self.label_2.setObjectName("label_2")
        self.label_2.setMinimumWidth(100)
        self.horizontalLayout.addWidget(self.label_2)
        self.textEdit = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget)
        self.textEdit.setObjectName("textEdit")
        self.textEdit.setDisabled(True)
        self.textEdit.setMaximumHeight(50)
        self.textEdit.textChanged.connect(self.checkIfFilled)
        self.horizontalLayout.addWidget(self.textEdit)
        self.pushButton = QtWidgets.QPushButton(parent=self.verticalLayoutWidget)
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.selectDeTable)
        self.horizontalLayout.addWidget(self.pushButton)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout2.setObjectName("horizontalLayout")
        self.label_3 = QtWidgets.QLabel(parent=self.verticalLayoutWidget)
        self.label_3.setObjectName("label_3")
        self.label_3.setMinimumWidth(100)
        self.horizontalLayout2.addWidget(self.label_3)
        self.textEdit_2 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget)
        self.textEdit_2.setObjectName("textEdit_2")
        self.textEdit_2.textChanged.connect(self.checkIfFilled)
        self.textEdit_2.setDisabled(True)
        self.textEdit_2.setMaximumHeight(50)
        self.horizontalLayout2.addWidget(self.textEdit_2)
        self.pushButton_3 = QtWidgets.QPushButton(parent=self.verticalLayoutWidget)
        self.pushButton_3.setObjectName("pushButton")
        self.pushButton_3.clicked.connect(self.browseDirectoryForOutput)
        self.horizontalLayout2.addWidget(self.pushButton_3)
        self.verticalLayout.addLayout(self.horizontalLayout2)
        self.pushButton_2 = QtWidgets.QPushButton(parent=self.verticalLayoutWidget)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.setMinimumHeight(30)
        self.pushButton_2.setDisabled(True)
        self.pushButton_2.clicked.connect(self.performGO)
        self.verticalLayout.addWidget(self.pushButton_2)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label.setText(_translate("Dialog", "GO Terms"))
        self.label_2.setText(_translate("Dialog", "DE Table"))
        self.label_3.setText(_translate("Dialog", "Output"))
        self.pushButton.setText(_translate("Dialog", "Browse"))
        self.pushButton_3.setText(_translate("Dialog", "Browse"))
        self.pushButton_2.setText(_translate("Dialog", "Start"))

if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    Group = QtWidgets.QMainWindow()
    ui = Ui_GO()
    ui.setupUi(Group)
    Group.show()
    sys.exit(app.exec())
