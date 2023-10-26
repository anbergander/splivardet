import subprocess
import time

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import QProcess

from cancel_run import Ui_CancelRun


class Ui_Progress(object):
    def runProcess(self):
        self.label.setText("Currently Running: Trimming")
        self.call()
        self.label.setText("Currently Running: Alignment")
        self.call2()
        self.label.setText("Currently Running: Isoform Detection")
        cmd = "nextflow run splivarquant"
        subprocess.call(cmd, shell=True)
    def call(self):
        cmd = "nextflow run porechop"
        subprocess.call(cmd, shell=True)
    def call2(self):
        cmd = "nextflow run splivardet"
        subprocess.call(cmd, shell=True)

    def openWindow(self):
        self.window = QtWidgets.QDialog()
        self.ui = Ui_CancelRun()
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
        self.pushButton.clicked.connect(self.openWindow)
        self.pushButton.setGeometry(QtCore.QRect(370, 80, 121, 41))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("cancel-vector-icon.jpg"), QtGui.QIcon.Mode.Normal,
                       QtGui.QIcon.State.Off)
        font.setPointSize(10)
        self.pushButton.setFont(font)
        self.pushButton.setIcon(icon)
        self.pushButton.setObjectName("pushButton")

        self.retranslateUi(Progress)
        QtCore.QMetaObject.connectSlotsByName(Progress)
        self.runProcess()

    def retranslateUi(self, Progress):
        _translate = QtCore.QCoreApplication.translate
        self.label.setText(_translate("runDialog", "<html><head/><body><p align=\"center\"><span style=\" font-weight:600;\">Currently running...</span></p></body></html>"))
        Progress.setWindowTitle(_translate("Progress", "Progress"))
        self.pushButton.setText(_translate("Progress", "Cancel"))


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    Progress = QtWidgets.QMainWindow()
    ui = Ui_Progress()
    ui.setupUi(Progress)
    Progress.show()
    sys.exit(app.exec())