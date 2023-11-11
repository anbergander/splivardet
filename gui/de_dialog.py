import os

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QFileDialog, QMessageBox

from de import de_final


class Ui_DE(object):
    def checkIfFilled(self):
        if (self.textEdit.toPlainText() == '' or self.textEdit_2.toPlainText() == '' or self.textEdit_3.toPlainText() == ''):
            self.pushButton_2.setDisabled(True)
        else:
            self.pushButton_2.setDisabled(False)
    def separateString(self):
        if self.textEdit_2.toPlainText() == '':
            return
        if not "," in self.textEdit_2.toPlainText():
            QMessageBox.information(self.verticalLayoutWidget, "Error", "Must be a List of Conditons e.g. C,V,C, ...")
        text = self.textEdit_2.toPlainText()
        li = list(text.split(","))
        return li
    def identifyContrasts(self, li):
        myset = set(li)
        if len(myset) > 2:
            QMessageBox.information(self.verticalLayoutWidget, "Error", "Currently DE is only available for two groups")
        return list(myset)
    def selectCountTable(self):
        fileName, _ = QFileDialog.getOpenFileName(Group, "Open CSV",
                                                  os.getenv('Home'), "CSV Files (*.csv);; TSV Files (*.tsv);; "
                                                                     "TXT Files (*.txt) ;; TAB Files (*.tab)")
        self.textEdit.setText(fileName)
    def browseDirectoryForOutput(self):
        fname = QFileDialog.getExistingDirectory(self.verticalLayoutWidget, "Select Directory", os.getenv('Home'))
        self.textEdit_3.setPlainText(fname)
    def performDe(self):
        li = self.separateString()
        contrasts = self.identifyContrasts(li)
        de_final.main(self.textEdit.toPlainText(), li, contrasts[0], contrasts[1], self.textEdit_3.toPlainText())
        QMessageBox.information(self.verticalLayoutWidget, "Done", "DE was performed")
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.setFixedSize(500, 256)
        self.label = QtWidgets.QLabel(parent=Dialog)
        self.label.setGeometry(QtCore.QRect(100, 10, 400, 20))
        self.label.setObjectName("label")
        font = QtGui.QFont()
        font.setPointSize(16)
        font.setFamily("Calibri")
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.verticalLayoutWidget = QtWidgets.QWidget(parent=Dialog)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(20, 50, 460, 183))
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
        self.textEdit.setMaximumHeight(50)
        self.textEdit.setDisabled(True)
        self.textEdit.textChanged.connect(self.checkIfFilled)
        self.horizontalLayout.addWidget(self.textEdit)
        self.pushButton = QtWidgets.QPushButton(parent=self.verticalLayoutWidget)
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.selectCountTable)
        self.horizontalLayout.addWidget(self.pushButton)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_3 = QtWidgets.QLabel(parent=self.verticalLayoutWidget)
        self.label_3.setObjectName("label_3")
        self.label_3.setMinimumWidth(100)
        self.horizontalLayout_2.addWidget(self.label_3)
        self.textEdit_2 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget)
        self.textEdit_2.setObjectName("textEdit_2")
        self.textEdit_2.setMaximumHeight(50)
        self.textEdit_2.textChanged.connect(self.checkIfFilled)
        self.textEdit_2.textChanged.connect(self.separateString)
        self.horizontalLayout_2.addWidget(self.textEdit_2)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout3.setObjectName("horizontalLayout3")
        self.label_4 = QtWidgets.QLabel(parent=self.verticalLayoutWidget)
        self.label_4.setObjectName("label_2")
        self.label_4.setMinimumWidth(100)
        self.horizontalLayout3.addWidget(self.label_4)
        self.textEdit_3 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget)
        self.textEdit_3.setObjectName("textEdit")
        self.textEdit_3.setMaximumHeight(50)
        self.textEdit_3.setDisabled(True)
        self.textEdit_3.textChanged.connect(self.checkIfFilled)
        self.horizontalLayout3.addWidget(self.textEdit_3)
        self.pushButton_3 = QtWidgets.QPushButton(parent=self.verticalLayoutWidget)
        self.pushButton_3.setObjectName("pushButton")
        self.pushButton_3.clicked.connect(self.browseDirectoryForOutput)
        self.horizontalLayout3.addWidget(self.pushButton_3)
        self.verticalLayout.addLayout(self.horizontalLayout3)
        self.pushButton_2 = QtWidgets.QPushButton(parent=self.verticalLayoutWidget)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.setMinimumHeight(30)
        self.pushButton_2.setDisabled(True)
        self.pushButton_2.clicked.connect(self.performDe)
        self.verticalLayout.addWidget(self.pushButton_2)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Start DEG"))
        self.label.setText(_translate("Dialog", "Differential Gene Expression"))
        self.label_2.setText(_translate("Dialog", "Count Table"))
        self.pushButton.setText(_translate("Dialog", "Browse"))
        self.pushButton_3.setText(_translate("Dialog", "Browse"))
        self.label_3.setText(_translate("Dialog", "Conditions"))
        self.label_3.setToolTip(_translate("Dialog", "Enter List of Conditions for each sample, separated by comma"))
        self.label_4.setText(_translate("Dialog", "Output"))
        self.textEdit_2.setPlaceholderText(_translate("Dialog", "Conditions for each sample, separated by comma"))
        self.pushButton_2.setText(_translate("Dialog", "Start"))

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Group = QtWidgets.QMainWindow()
    ui = Ui_DE()
    ui.setupUi(Group)
    Group.show()
    sys.exit(app.exec())