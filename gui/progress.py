import subprocess

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QMessageBox

from cancel_run import Ui_CancelRun
from gui.report import Ui_Report


class Ui_Progress(object):
    def runProcess(self):
        self.pushButton_2.hide()
        self.call()
        self.call2()
        self.call3()
        self.pushButton.show()
        QMessageBox.information(self.window, "Done", "Isoform Detection was performed")
    def call(self):
        cmd = "nextflow run porechop"
        subprocess.call(cmd, shell=True)
    def call2(self):
        cmd = "nextflow run splivardet"
        subprocess.call(cmd, shell=True)

    def call3(self):
        cmd = "nextflow run splivarquant"
        subprocess.call(cmd, shell=True)

    def openWindow(self):
        self.window = QtWidgets.QDialog()
        self.ui = Ui_CancelRun()
        self.ui.setupUi(self.window)
        self.window.show()

    def openReport(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_Report()
        self.ui.setupUi(self.window)
        self.window.show()

    def setupUi(self, Progress):
        self.window = self.window = QtWidgets.QMainWindow()
        Progress.setObjectName("Progress")
        Progress.setFixedSize(500, 155)
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setBold(True)
        font.setWeight(75)
        Progress.setFont(font)

        self.label = QtWidgets.QLabel(parent=Progress)
        self.label.setGeometry(QtCore.QRect(70, 20, 400, 31))
        self.label.setFont(font)
        self.label.setObjectName("label")

        self.pushButton = QtWidgets.QPushButton(parent=Progress)
        self.pushButton.clicked.connect(self.openReport)
        self.pushButton.clicked.connect(Progress.close)
        self.pushButton.setGeometry(QtCore.QRect(370, 80, 121, 41))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("next-step-28.png"), QtGui.QIcon.Mode.Normal,
                       QtGui.QIcon.State.Off)
        font.setPointSize(10)
        self.pushButton.setFont(font)
        self.pushButton.setIcon(icon)
        self.pushButton.setObjectName("pushButton")
        self.pushButton.hide()
        self.pushButton_2 = QtWidgets.QPushButton(parent=Progress)
        self.pushButton_2.clicked.connect(self.pushButton_2.hide)
        self.pushButton_2.clicked.connect(self.runProcess)
        self.pushButton_2.setGeometry(QtCore.QRect(15, 80, 121, 41))
        self.pushButton_2.setFont(font)

        self.retranslateUi(Progress)
        QtCore.QMetaObject.connectSlotsByName(Progress)

    def retranslateUi(self, Progress):
        _translate = QtCore.QCoreApplication.translate
        self.label.setText(_translate("runDialog", "<html><head/><body><p align=\"center\"><span style=\" font-weight:600;\">Currently running...</span></p></body></html>"))
        Progress.setWindowTitle(_translate("Progress", "Progress"))
        self.pushButton.setText(_translate("Progress", "Continue"))
        self.pushButton_2.setText(_translate("Progress", "Start"))


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    Progress = QtWidgets.QMainWindow()
    ui = Ui_Progress()
    ui.setupUi(Progress)
    Progress.show()
    sys.exit(app.exec())