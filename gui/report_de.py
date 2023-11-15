import os
import sys
import textwrap
from datetime import datetime

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import pyqtSignal, QObject, QRunnable, QThreadPool, pyqtSlot
from PyQt6.QtWidgets import QFileDialog, QMessageBox
from pdfrw import PdfReader
from pdfrw.buildxobj import pagexobj
from pdfrw.toreportlab import makerl
from reportlab.pdfgen.canvas import Canvas

class WorkerSignals(QObject):
    error = pyqtSignal(str)
    file_saved_as = pyqtSignal(str)

class PdfGenerator(QRunnable):
    def __init__(self, data, folder):
        super().__init__()
        self.folder = folder
        self.data = data
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        try:
            outfile = self.folder + "//" + "result_de.pdf"

            template = PdfReader("report_de.pdf", decompress=False).pages[0]
            template_obj = pagexobj(template)

            canvas = Canvas(outfile)

            xobj_name = makerl(canvas, template_obj)
            canvas.doForm(xobj_name)

            ystart = 690

            today = datetime.today()
            canvas.drawString(450, ystart, today.strftime('%F'))

            canvas.drawString(180, ystart - 65, self.data['experiment'])
            canvas.drawString(180, ystart - 89, self.data['reference'])
            canvas.drawString(180, ystart - 113, self.data['gtf'])


            canvas.drawString(235, ystart - 190, self.data['de_method'])
            de_output = self.data['de_out'].split('\n')
            n=0
            if de_output:
                for elem in de_output:
                    canvas.drawString(235, 460 - (n * 28), elem)
                    n= n+1

            canvas.drawString(235, ystart - 450, self.data['gsea_methode'])
            gsea_output = self.data['gsea_out'].split('\n')
            n = 0
            if gsea_output:
                for elem in gsea_output:
                    canvas.drawString(235, 200 - (n * 28), elem)
                    n = n + 1

            canvas.save()

        except Exception as e:
            self.signals.error.emit(str(e))
            return

        self.signals.file_saved_as.emit(outfile)

class Ui_ReportDe(object):
    def generate(self):
        fname = QFileDialog.getExistingDirectory(self.centralwidget, "Select Directory", os.getenv('Home'))
        fname = fname.replace('\\', '/')
        self.path = fname
        if fname:
            self.pushButton.setDisabled(True)
            data = {
                'experiment': self.textEdit_3.toPlainText(),
                'gtf': self.textEdit_5.toPlainText(),
                'reference': self.textEdit_4.toPlainText(),
                'de_method': self.textEdit.toPlainText(),
                'de_out': self.textBrowser.toPlainText(),
                'gsea_methode': self.textEdit_2.toPlainText(),
                'gsea_out': self.textBrowser_2.toPlainText(),
            }
            g = PdfGenerator(data, fname)
            g.signals.file_saved_as.connect(self.generated)
            g.signals.error.connect(print)
            self.threadpool.start(g)
    def generated(self, outfile):
        self.pushButton.setDisabled(True)
        try:
            os.startfile(outfile)
        except Exception:
            QMessageBox.information(self.centralwidget, "Finished", "PDF has been generated")

    def setupUi(self, MainWindow):
        self.path = ""
        self.threadpool = QThreadPool()
        MainWindow.setObjectName("MainWindow")
        MainWindow.setFixedSize(875, 370)
        self.centralwidget = QtWidgets.QWidget(parent=MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.label = QtWidgets.QLabel(parent=self.centralwidget)
        self.label.setGeometry(QtCore.QRect(220, 0, 420, 31))
        self.label.setObjectName("label")
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        self.label.setFont(font)
        self.verticalLayoutWidget_3 = QtWidgets.QWidget(parent=self.centralwidget)
        self.verticalLayoutWidget_3.setGeometry(QtCore.QRect(20, 40, 830, 300))
        self.verticalLayoutWidget_3.setObjectName("verticalLayoutWidget_3")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_3)
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.label_8 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_3)
        self.label_8.setObjectName("label_8")
        self.label_8.setMinimumWidth(100)
        self.horizontalLayout_7.addWidget(self.label_8)
        self.textEdit_3 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_3)
        self.textEdit_3.setObjectName("textEdit_3")
        self.textEdit_3.setMaximumHeight(30)
        self.horizontalLayout_7.addWidget(self.textEdit_3)
        self.verticalLayout_3.addLayout(self.horizontalLayout_7)
        self.horizontalLayout_11 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.label_9 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_3)
        self.label_9.setObjectName("label_9")
        self.label_9.setMinimumWidth(100)
        self.horizontalLayout_8.addWidget(self.label_9)
        self.textEdit_4 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_3)
        self.textEdit_4.setObjectName("textEdit_4")
        self.textEdit_4.setMaximumHeight(30)
        self.horizontalLayout_8.addWidget(self.textEdit_4)
        self.horizontalLayout_11.addLayout(self.horizontalLayout_8)
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.label_10 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_3)
        self.label_10.setObjectName("label_10")
        self.label_10.setMinimumWidth(100)
        self.horizontalLayout_9.addWidget(self.label_10)
        self.textEdit_5 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_3)
        self.textEdit_5.setObjectName("textEdit_5")
        self.textEdit_5.setMaximumHeight(30)
        self.horizontalLayout_9.addWidget(self.textEdit_5)
        self.horizontalLayout_11.addLayout(self.horizontalLayout_9)
        self.verticalLayout_3.addLayout(self.horizontalLayout_11)
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.label_2 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_3)
        self.label_2.setObjectName("label_2")
        self.label_2.setMinimumWidth(100)
        self.horizontalLayout_3.addWidget(self.label_2)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_3 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_3)
        self.label_3.setObjectName("label_3")
        self.label_3.setMinimumWidth(100)
        self.horizontalLayout.addWidget(self.label_3)
        self.textEdit = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_3)
        self.textEdit.setObjectName("textEdit")
        self.textEdit.setMaximumHeight(30)
        self.textEdit.setDisabled(True)
        self.horizontalLayout.addWidget(self.textEdit)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_4 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_3)
        self.label_4.setObjectName("label_4")
        self.label_4.setMinimumWidth(100)
        self.horizontalLayout_2.addWidget(self.label_4)
        self.textBrowser = QtWidgets.QTextBrowser(parent=self.verticalLayoutWidget_3)
        self.textBrowser.setObjectName("textBrowser")
        self.textBrowser.setMaximumHeight(150)
        self.textBrowser.setDisabled(True)
        self.horizontalLayout_2.addWidget(self.textBrowser)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_10.addLayout(self.verticalLayout)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem2)
        self.label_5 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_3)
        self.label_5.setObjectName("label_5")
        self.label_5.setMinimumWidth(100)
        self.horizontalLayout_4.addWidget(self.label_5)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem3)
        self.verticalLayout_2.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_6 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_3)
        self.label_6.setObjectName("label_6")
        self.label_6.setMinimumWidth(100)
        self.horizontalLayout_5.addWidget(self.label_6)
        self.textEdit_2 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_3)
        self.textEdit_2.setObjectName("textEdit_2")
        self.textEdit_2.setMaximumHeight(30)
        self.textEdit_2.setDisabled(True)
        self.horizontalLayout_5.addWidget(self.textEdit_2)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_7 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_3)
        self.label_7.setObjectName("label_7")
        self.label_7.setMinimumWidth(100)
        self.horizontalLayout_6.addWidget(self.label_7)
        self.textBrowser_2 = QtWidgets.QTextBrowser(parent=self.verticalLayoutWidget_3)
        self.textBrowser_2.setObjectName("textBrowser_2")
        self.textBrowser_2.setMaximumHeight(150)
        self.textBrowser_2.setDisabled(True)
        self.horizontalLayout_6.addWidget(self.textBrowser_2)
        self.verticalLayout_2.addLayout(self.horizontalLayout_6)
        self.horizontalLayout_10.addLayout(self.verticalLayout_2)
        self.verticalLayout_3.addLayout(self.horizontalLayout_10)
        self.pushButton = QtWidgets.QPushButton(parent=self.verticalLayoutWidget_3)
        self.pushButton.setObjectName("pushButton")
        self.pushButton.setMinimumHeight(30)
        self.pushButton.clicked.connect(self.generate)
        self.verticalLayout_3.addWidget(self.pushButton)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(parent=MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 571, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(parent=MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Report"))
        self.label.setText(_translate("MainWindow", "Report: DEG and Pathway Enrichment"))
        self.label_8.setText(_translate("MainWindow", "Experiment"))
        self.textEdit_3.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.label_9.setText(_translate("MainWindow", "Reference genome"))
        self.textEdit_4.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.label_10.setText(_translate("MainWindow", "Annotation"))
        self.textEdit_5.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.label_2.setText(_translate("MainWindow", "DEG"))
        self.label_3.setText(_translate("MainWindow", "Method"))
        self.textEdit.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">PyDESeq2</span></p></body></html>"))
        self.label_4.setText(_translate("MainWindow", "Output"))
        self.textBrowser.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">Metadata (metadata.tsv)</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\"> DE Table (result_table.tsv)</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\"> Significant Genes (significant_de.tsv)</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\"> Volcano-Plot (volcano.svg; volcano.png)</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\"> Heatmap (heatmap_de.svg)</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\"> PCA (pca_normalized.svg)</span></p></body></html>"))
        self.label_5.setText(_translate("MainWindow", "Pathway Enrichment"))
        self.label_6.setText(_translate("MainWindow", "Method"))
        self.textEdit_2.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\">GSEApy</span></p></body></html>"))
        self.label_7.setText(_translate("MainWindow", "Output"))
        self.textBrowser_2.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:\'Segoe UI\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\"> GSEA Table (gsea.tsv)</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\"> Enrichment Plots (gsea_enrichment_plot.svg)</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\"> Enriched Pathways (gsea_terms_enriched.svg)</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt;\"> Reduced Pathways (gsea_terms_decreased.svg)</span></p></body></html>"))
        self.pushButton.setText(_translate("MainWindow", "Generate"))

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    Group = QtWidgets.QMainWindow()
    ui = Ui_ReportDe()
    ui.setupUi(Group)
    Group.show()
    sys.exit(app.exec())
