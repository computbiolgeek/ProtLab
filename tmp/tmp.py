

class AudioFile:
    """
    """
    def __init__( self, filename ):
        if not filename.lower().endswith( self.ext ):
            raise Exception( "Invalid file format" )
        self.filename = filename

class MP3File( AudioFile ):
    """
    """
    ext = "mp3"
    def play( self ):
        print( "playing {} as mp3".format( self.filename ) )

class WavFile( AudioFile ):
    """
    """
    ext = "wav"
    def play( self ):
        print( "playing {} as wav".format( self.filename ) )

class OggFile( AudioFile ):
    """
    """
    ext = "ogg"
    def play( self ):
        print( "playing {} as ogg".format( self.filename ) )

def PlayAudio( audio_file ):
    audio_file.play()

class FlacFile:
    """
    """
    def __init__( self, filename ):
        if not filename.lower().endswith( "flac" ):
            raise Exception( "Invalid file format" )
        self.filename = filename

    def play( self ):
        print( "playing {} as flac".format( self.filename ) )


class InvalidWithdrawal(Exception):
    """

    """
    def __init__(self, balance, amount):
        super().__init__("Account does not have ${}".format(amount))
        self._balance = balance
        self._amount = amount

    def overage(self):
        return self._amount - self._balance

try:
    raise InvalidWithdrawal(25, 50)
except InvalidWithdrawal as e:
    print("I'm sorry, but your withdrawal is more than your balance "
          "by ${}".format(e.overage()))


class Silly:
    @property
    def silly(self):
        "This is is silly property"
        print("You're getting silly")
        return self._silly

    @silly.setter
    def silly(self, value):
        print("You're are making silly {}".format(value))
        self._silly = value

    @silly.deleter
    def silly(self):
        print("Whoah, you killed silly!")
        del self._silly


from urllib.request import urlopen

class WebPage:
    def __init__(self, url):
        self.url = url
        self._content = None

    @property
    def content(self):
        if not self._content:
            print("Retrieving New Page ...")
            self._content = urlopen(self.url).read()
        return self._content