# an object of WSGI application
import os

from flask import Flask, render_template
print(os.getcwd())
app = Flask(__name__, template_folder='./templates')  # Flask constructor


# A decorator used to tell the application
# which URL is associated function
@app.route('/')
def hello():
    my_dir = os.getcwd()

    return render_template("index.html")


if __name__ == '__main__':
    app.debug = True
    app.run()
