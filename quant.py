from scipy.signal import savgol_filter


class Quant_smooth:
    def __init__(self, window):
        self.BOTTOM_REL_INTEN = 0.1
        self.MAX_REACH_BOTTOM_TIME = 2
        self.MAX_UP_OR_DOWN = 2
        self.LOCAL_LOOK_WINDOW = 5
        self.smooth_window = window

    def find_maxima_left(self, intens, gpsm_idx):
        i = gpsm_idx
        max_idx = gpsm_idx
        down_num = 0
        global_min = 1e100
        local_min = 1e100
        local_look_num = 0
        while i >= 0:
            if intens[i] <= 0:
                return max_idx
            elif intens[i] > intens[max_idx]:
                max_idx = i
            if intens[i] < local_min:
                local_min = intens[i]
            local_look_num += 1
            if local_look_num == self.LOCAL_LOOK_WINDOW:
                if local_min < global_min:
                    global_min = local_min
                    down_num += 1
                    if down_num > self.MAX_UP_OR_DOWN:
                        return max_idx
                local_min = 1e100
                local_look_num = 0
            i -= 1
        return max_idx

    def find_maxima_right(self, intens, gpsm_idx):
        i = gpsm_idx
        max_idx = gpsm_idx
        down_num = 0
        global_min = 1e100
        local_min = 1e100
        local_look_num = 0
        while i < len(intens):
            if intens[i] <= 0:
                return max_idx
            elif intens[i] > intens[max_idx]:
                max_idx = i
            if intens[i] < local_min:
                local_min = intens[i]
            local_look_num += 1
            if local_look_num == self.LOCAL_LOOK_WINDOW:
                if local_min < global_min:
                    global_min = local_min
                    down_num += 1
                    if down_num > self.MAX_UP_OR_DOWN:
                        return max_idx
                local_min = 1e100
                local_look_num = 0
            i += 1
        return max_idx

    def local_maxima_pos(self, intens, gpsm_idx):
        left_max_idx = self.find_maxima_left(intens, gpsm_idx)
        right_max_idx = self.find_maxima_right(intens, gpsm_idx)
        return (
            left_max_idx
            if intens[left_max_idx] > intens[right_max_idx]
            else right_max_idx
        )

    def Retention_left(self, intens, gpsm_idx, max_pos):
        bottom_inten = intens[max_pos] * self.BOTTOM_REL_INTEN
        bottom_time = 0
        left = (gpsm_idx - 1) if max_pos > gpsm_idx else (max_pos - 1)
        up_num = 0
        global_max = 0
        local_max = 0
        local_look_num = 0
        while left >= 0:
            if intens[left] <= 0 or intens[left] > intens[max_pos]:
                break
            elif intens[left] < bottom_inten:
                bottom_time += 1
                if bottom_time == self.MAX_REACH_BOTTOM_TIME:
                    break
            if intens[left] > local_max:
                local_max = intens[left]
            local_look_num += 1
            if local_look_num == self.LOCAL_LOOK_WINDOW:
                if local_max > global_max:
                    global_max = local_max
                    up_num += 1
                    if up_num > self.MAX_UP_OR_DOWN:
                        break
                local_max = 0
                local_look_num = 0
            left -= 1
        left = left if left >= 0 else 0
        return left

    def Retention_right(self, intens, gpsm_idx, max_pos):
        bottom_inten = intens[max_pos] * self.BOTTOM_REL_INTEN
        bottom_num = 0
        right = (gpsm_idx + 1) if max_pos < gpsm_idx else (max_pos + 1)
        up_num = 0
        global_max = 0
        local_max = 0
        local_look_num = 0
        while right < len(intens):
            if intens[right] <= 0 or intens[right] > intens[max_pos]:
                break
            elif intens[right] < bottom_inten:
                bottom_num += 1
                if bottom_num == self.MAX_REACH_BOTTOM_TIME:
                    break
            if intens[right] > local_max:
                local_max = intens[right]
            local_look_num += 1
            if local_look_num == self.LOCAL_LOOK_WINDOW:
                if local_max > global_max:
                    global_max = local_max
                    up_num += 1
                    if up_num > self.MAX_UP_OR_DOWN:
                        break
                local_max = 0
                local_look_num = 0
            right += 1
        right = right if right < len(intens) else (len(intens) - 1)
        return right

    def Retention_range_smooth(self, intens, gpsm_idx, smooth_window):
        smooth_intens = savgol_filter(
            intens, window_length=smooth_window, polyorder=2, mode="constant", cval=0
        )
        smooth_max_pos = self.local_maxima_pos(smooth_intens, gpsm_idx)
        left_start = gpsm_idx - 1
        left = self.Retention_left(smooth_intens, left_start, smooth_max_pos)
        right_start = gpsm_idx + 1
        right = self.Retention_right(smooth_intens, right_start, smooth_max_pos)

        real_max_pos = smooth_max_pos
        i = smooth_max_pos - self.LOCAL_LOOK_WINDOW
        end = smooth_max_pos + self.LOCAL_LOOK_WINDOW + 1
        if i < 0:
            i = 0
        if end >= len(intens):
            end = len(intens)
        for i in range(i, end):
            if intens[i] > intens[real_max_pos]:
                real_max_pos = i

        if left > real_max_pos:
            left = real_max_pos - 1
        elif left < 0:
            left = 0
        while left >= 1 and intens[left - 1] > 0:
            if intens[left - 1] > intens[left]:
                break
            left -= 1
        if right < real_max_pos:
            right = real_max_pos + 1
        elif right >= len(intens):
            right = len(intens) - 1
        while right < len(intens) - 1 and intens[right + 1] > 0:
            if intens[right + 1] > intens[right]:
                break
            right += 1
        return left, right, real_max_pos

    def side_trim(self, intens, gpsm_idx, left, right, max_pos):
        left_bound = gpsm_idx if gpsm_idx < max_pos else max_pos
        right_bound = gpsm_idx if gpsm_idx > max_pos else max_pos
        left_min = 1e100
        right_min = 1e100
        if gpsm_idx >= 1 and intens[gpsm_idx - 1] > 0 and left_bound == gpsm_idx:
            left_bound = gpsm_idx - 1
        if (
            gpsm_idx + 1 < len(intens)
            and intens[gpsm_idx + 1] > 0
            and right_bound == gpsm_idx
        ):
            right_bound = gpsm_idx + 1
        for i in range(left_bound, left - 1, -1):
            if intens[i] == 0:
                left = i + 1
                break
        for i in range(right_bound, right + 1):
            if intens[i] == 0:
                right = i - 1
                break
        if left > gpsm_idx:
            left = gpsm_idx
        if right < gpsm_idx:
            right = gpsm_idx

        # for i in range(left, left_bound+1):
        #    if intens[i] > 0 and intens[i] < left_min: left_min = intens[i]
        # for i in range(right_bound, right+1):
        #    if intens[i] > 0 and intens[i] < right_min: right_min = intens[i]
        # while left < left_bound and intens[left]>left_min+1e-6: left += 1
        # while right > right_bound and intens[right]>right_min+1e-6: right -= 1

        return left, right, max_pos

    def Retention_range(self, intens, gpsm_idx):
        left, right, max_pos = self.Retention_range_smooth(
            intens, gpsm_idx, self.smooth_window
        )

        return self.side_trim(intens, gpsm_idx, left, right, max_pos)
