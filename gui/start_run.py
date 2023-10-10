from PyQt6 import QtCore, QtGui, QtWidgets

from progress import Ui_Progress


class Ui_runDialog(object):
    def openWindow(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_Progress()
        self.ui.setupUi(self.window)
        self.window.show()
    def setupUi(self, runDialog):
        runDialog.setObjectName("runDialog")
        runDialog.setFixedSize(470, 150)
        self.label = QtWidgets.QLabel(parent=runDialog)
        self.label.setGeometry(QtCore.QRect(40, 0, 400, 61))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.pushButton_3 = QtWidgets.QPushButton(parent=runDialog)
        self.pushButton_3.setGeometry(QtCore.QRect(20, 70, 121, 41))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.pushButton_3.setFont(font)
        self.pushButton_3.clicked.connect(runDialog.close)
        self.pushButton = QtWidgets.QPushButton(parent=runDialog)
        self.pushButton.setGeometry(QtCore.QRect(330, 70, 121, 41))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(10)
        self.pushButton.setFont(font)
        self.pushButton.clicked.connect(runDialog.close)
        self.pushButton.clicked.connect(self.openWindow)
        self.pushButton.setObjectName("pushButton")

        self.retranslateUi(runDialog)
        QtCore.QMetaObject.connectSlotsByName(runDialog)

    def retranslateUi(self, runDialog):
        _translate = QtCore.QCoreApplication.translate
        runDialog.setWindowTitle(_translate("runDialog", "Start Run"))
        self.label.setText(_translate("runDialog", "<html><head/><body><p align=\"center\"><span style=\" font-weight:600;\">Do you want to start your run?</span></p></body></html>"))
        self.pushButton_3.setText(_translate("runDialog", "No"))
        self.pushButton.setText(_translate("runDialog", "Yes"))

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Window = QtWidgets.QMainWindow()
    ui = Ui_runDialog()
    ui.setupUi(Window)
    Window.show()
    sys.exit(app.exec())