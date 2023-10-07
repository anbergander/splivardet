from PyQt6 import QtCore, QtGui, QtWidgets


class Ui_CancelRun(object):

    def setupUi(self, CancelRun):
        self.window = self.window = QtWidgets.QMainWindow()
        CancelRun.setObjectName("ParameterFormReport")
        CancelRun.resize(400, 200)
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setBold(True)
        font.setWeight(75)
        CancelRun.setFont(font)

        self.buttonBox = QtWidgets.QDialogButtonBox(parent=CancelRun)
        self.buttonBox.setGeometry(QtCore.QRect(40, 150, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Orientation.Horizontal)
        self.buttonBox.setStandardButtons(
            QtWidgets.QDialogButtonBox.StandardButton.No | QtWidgets.QDialogButtonBox.StandardButton.Yes)
        self.buttonBox.setObjectName("buttonBox")
        self.textBrowser = QtWidgets.QTextBrowser(parent=CancelRun)
        self.textBrowser.setGeometry(QtCore.QRect(60, 10, 291, 121))
        self.textBrowser.setObjectName("textBrowser")

        self.retranslateUi(CancelRun)
        self.buttonBox.accepted.connect(CancelRun.accept)  # type: ignore
        self.buttonBox.rejected.connect(CancelRun.reject)  # type: ignore
        QtCore.QMetaObject.connectSlotsByName(CancelRun)

        self.retranslateUi(CancelRun)
        QtCore.QMetaObject.connectSlotsByName(CancelRun)

    def retranslateUi(self, CancelRun):
        _translate = QtCore.QCoreApplication.translate
        CancelRun.setWindowTitle(_translate("CancelRun", "Cancel Run"))
        self.textBrowser.setHtml(_translate("CancelRun",
                                            "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
                                            "<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
                                            "p, li { white-space: pre-wrap; }\n"
                                            "</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
                                            "<p align=\"center\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:16pt;\">Are you sure you want to cancel your run?</span></p>\n"
                                            "<p align=\"center\" style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:16pt;\"><br /></p>\n"
                                            "<p align=\"center\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:16pt;\">Data might be lost.</span></p></body></html>"))


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    CancelRun = QtWidgets.QDialog()
    ui = Ui_CancelRun()
    ui.setupUi(CancelRun)
    CancelRun.show()
    sys.exit(app.exec())