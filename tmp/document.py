"""

"""

class Document:
    """

    """
    def __init__(self):
        self.filename = ''
        self.title = ''
        self.characters = []
        self.cursor = Cursor(self)

    def insert(self, character):
        self.characters.insert(self.cursor, character)
        self.cursor += 1

    def delete(self):
        del self.characters[self.cursor]

    def clear(self):
        self.characters.clear()
        self.cursor = 0

    def save(self):
        with open(self.filename, "wt") as f:
            f.write("".join(self.characters))

    def forward(self):
        self.cursor += 1

    def back(self):
        self.cursor -= 1


class Cursor:
    """

    """
    def __init__(self, document):
        self.document = document
        self.position = 0

    def forward(self):
        self.position += 1

    def back(self):
        self.position -= 1

    def home(self):
        if self.position == 0:
            return
        while self.document.characters[self.position - 1] != '\n':
            self.back()
            if self.position == 0:
                break

    def end(self):
        while self.position < len(self.document.characters)
            and self.document.characters[self.position - 1] != '\n':
            self.forward()