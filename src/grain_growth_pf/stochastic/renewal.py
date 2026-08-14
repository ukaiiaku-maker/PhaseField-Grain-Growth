from dataclasses import dataclass


@dataclass
class RenewalWindow:
    length: float
    progress: float = 0.0

    def advance(self, geometric_increment: float) -> bool:
        self.progress += abs(geometric_increment)
        if self.progress >= self.length:
            self.progress %= self.length
            return True
        return False

