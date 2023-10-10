from PyQt6 import QtCore, QtGui, QtWidgets
from cancel_run import Ui_CancelRun


class Ui_Progress(object):
    def openWindow(self):
        self.window = QtWidgets.QDialog()
        self.ui = Ui_CancelRun()
        self.ui.setupUi(self.window)
        self.window.show()

    def setupUi(self, Progress):
        self.window = self.window = QtWidgets.QMainWindow()
        Progress.setObjectName("Progress")
        Progress.setFixedSize(400, 155)
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setBold(True)
        font.setWeight(75)
        Progress.setFont(font)

        self.label = QtWidgets.QLabel(parent=Progress)
        self.label.setGeometry(QtCore.QRect(70, 20, 261, 31))
        self.label.setFont(font)
        self.label.setObjectName("label")

        self.pushButton = QtWidgets.QPushButton(parent=Progress)
        self.pushButton.clicked.connect(self.openWindow)
        self.pushButton.setGeometry(QtCore.QRect(260, 80, 121, 41))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("cancel-vector-icon.jpg"), QtGui.QIcon.Mode.Normal,
                       QtGui.QIcon.State.Off)
        font.setPointSize(10)
        self.pushButton.setFont(font)
        self.pushButton.setIcon(icon)
        self.pushButton.setObjectName("pushButton")

        self.retranslateUi(Progress)
        QtCore.QMetaObject.connectSlotsByName(Progress)

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