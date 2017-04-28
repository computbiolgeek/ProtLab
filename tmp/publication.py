"""

"""


class Publication(object):
    """

    """
    def __init__(self, title, publisher, year):
        self._title = title
        self._publisher = publisher
        self._year = year

    def get_title(self):
        return self._title

    def get_publisher(self):
        return self._publisher

    def get_year(self):
        return self._year

    def set_title(self, title):
        self._title = title

    def set_publisher(self, publisher):
        self._publisher = publisher

    def set_year(self, year):
        self._year = year

    def display_info(self):
        print("Title: {}".format(self._title))
        print("Publisher: {}".format(self._publisher))
        print("Published in: {}".format(self._year))


class Book(Publication):
    """

    """
    def __init__(self, title, publisher, author, year, edition, isbn, price):
        super().__init__(title, publisher, year)
        self._author = author
        self._edition = edition
        self._isbn = isbn
        self._price = price

    def get_author(self):
        return self._author

    def get_edition(self):
        return self._edition

    def get_isbn(self):
        return self._isbn

    def get_price(self):
        return self._price

    def set_author(self, author):
        self._author = author

    def set_edition(self, edition):
        self._edition = edition

    def set_isbn(self, isbn):
        self._isbn = isbn

    def set_price(self, price):
        self._price = price

    def display_info(self):
        super().display_info()
        print("Author: {}".format(self._author))
        print("Edition: {}".format(self._edition))
        print("ISBN: {}".format(self._isbn))
        print("Price: {}".format(self._price))


class Magazine(Publication):
    """

    """
    def __init__(self, title, publisher, year, type, issue, rate):
        super().__init__(title, publisher, year)
        self._type = type
        self._issue = issue
        self._rate = rate