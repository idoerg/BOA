
from nltk import PorterStemmer
"""Format the text to aid naive bayes"""
def formatText(text):
    text = text.lower()
    text = text.replace('.',' ')
    text = text.replace('\\',' ')
    text = text.replace('/',' ')
    text = text.replace('\"',' ')
    text = text.replace('\'',' ')
    text = text.replace(':',' ')
    text = text.replace(';',' ')
    text = text.replace('(',' ')
    text = text.replace(')',' ')
    porter = PorterStemmer()
    return ' '.join([porter.stem(word) for word in text.split(' ')])

