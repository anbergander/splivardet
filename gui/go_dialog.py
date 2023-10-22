import os

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QFileDialog


class Ui_GO(object):
    def selectDeTable(self):
        fileName, _ = QFileDialog.getOpenFileName(self.verticalLayoutWidget, "Open CSV",
                                                  os.getenv('Home'), "CSV Files (*.csv);; TSV Files (*.tsv);; "
                                                                     "TXT Files (*.txt) ;; TAB Files (*.tab)")
        self.textEdit.setPlaceholderText(fileName)
        self.pushButton_2.setDisabled(False)

    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(470, 175)
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
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(20, 20, 430, 150))
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
        self.horizontalLayout.addWidget(self.textEdit)
        self.pushButton = QtWidgets.QPushButton(parent=self.verticalLayoutWidget)
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.selectDeTable)
        self.horizontalLayout.addWidget(self.pushButton)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.pushButton_2 = QtWidgets.QPushButton(parent=self.verticalLayoutWidget)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.setMinimumHeight(30)
        self.pushButton_2.setDisabled(True)
        self.verticalLayout.addWidget(self.pushButton_2)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label.setText(_translate("Dialog", "GO Terms"))
        self.label_2.setText(_translate("Dialog", "DE Table"))
        self.pushButton.setText(_translate("Dialog", "Browse"))
        self.pushButton_2.setText(_translate("Dialog", "Start"))

if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    Group = QtWidgets.QMainWindow()
    ui = Ui_GO()
    ui.setupUi(Group)
    Group.show()
    sys.exit(app.exec())
