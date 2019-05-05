import utils
import bwt_build

class FMSimpleIndex(object):
    def __init__(self, data, bw):
        self.data = bw.transform(data)
        self.offset = {}
        self._build(data)

    def _build(self, data):
        """ build the index """
        self.occ = utils.calc_first_occ(self.data)

    def _occ(self, qc):
        """ get the first occurance of letter qc in left-column"""
        c = self.occ.get(qc)
        if c == None:
            return 0
        return c

    def _count(self, idx, qc):
        """ count the occurances of letter qc (rank of qc) upto position idx """
        if not qc in self.occ.keys(): return 0
        c = 0
        for i in range(idx):
            if self.data[i] == qc:
                c += 1
        return c

    def _lf(self, idx, qc):
        """ get the nearset lf mapping for letter qc at position idx """
        o = self._occ(qc)
        c = self._count(idx, qc)
        return o + c

    def _walk(self, idx):
        """ find the offset in position idx of transformed string
            from the beginning """

        # walk to the beginning using lf mapping
        # this is same as inverse of burrow wheeler transformation
        # from arbitrary location
        r = 0
        i = idx
        while self.data[i] != self.bw.EOS:
            if self.offset.get(i):
                # we have cached the location and can use it
                r += self.offset[i]
                break
            r += 1
            i = self._lf(i, self.data[i])

        # save the offset of some idx for faster searches
        if not self.offset.get(idx):
            self.offset[i] = r
        return r

    def suffix(self, i):
        count = 0
        while self.sa[i] == None:
            i = self._lf(i, self.data[i])
            count = count + 1
        return self.sa[i] + count

    def bounds(self, q):
        """ find the first and last suffix positions for query q """
        top = 0
        bot = len(self.data)
        for i, qc in enumerate(q[::-1]):
            top = self._lf(top, qc)
            bot = self._lf(bot, qc)
            if top == bot: return (-1, -1)
        return (top, bot)

    def search(self, q):
        """ search the positions of query q """

        # find the suffixes for the query
        top, bot = self.bounds(q)
        matches = []
        # find the location of the suffixes
        # by walking the reverse text from that position
        # with lf mapping
        for i in range(top, bot):
            # pos = self._walk(i)
            pos = self.suffix(i)
            matches.append(pos)
        return sorted(matches)

    def count(self, q):
        """ count occurances of q in the index """
        top, bot = self.bounds(q)
        return bot - top

class FMCheckpointing(FMSimpleIndex):
    """ creates LF index with checkpoints """

    def __init__(self, data, bw, sa_step = 1, tally_step = 1):
        bwt_ref, idx, self.sa = bw.bwt(data)
        encoded = "".join(data[x] for x in bwt_ref)
        self.data = encoded
        self.offset = {}
        self.step = tally_step
        self.step_sa = sa_step
        self._build()

    def _build(self):
        """ build the index """
        self.occ = utils.calc_first_occ(self.data)
        self.C = utils.calc_checkpoints(self.data, self.step)
        self.sa = utils.calc_sa_checkpoints(self.sa, self.step_sa)

    def _count(self, idx, qc):
        """ count the occurances of letter qc (rank of qc) upto position idx """
        count = utils.count_letter_with_checkpoints(self.C, self.step, self.data, idx, qc)
        return count