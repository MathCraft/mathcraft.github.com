[latexpage]大家在学习动态规划的时候也许都接触过求[最长公共子序列(LCS)][1]算法，如果把原子类型的元素的相等性测试的代价看做是$O(1)$的话， LCS 经典的动态规划算法的时间复杂度为 $O(n^2)$， 这里的 $n$ 表示的是两个被对比序列的最大长度。

大家如果对[最长公共子序列][1]感兴趣的话，互联网上有许多相关的资料， 至少在这篇文章里我们不会深入这个话题。 我们的今天主题和这个算法有很密切的关系， 因为我们将使用 Mathematica 内置的函数 [`SequenceAlignment`][2] 来构建一个 Mathematica 版本的 [`diff`][3] 函数， 而 [`SequenceAlignment`][2] 背后使用的算法就恰恰是 LCS 的经典动态规划算法 (事实上也有 $O(n \log n)$ 的算法，在后来的版本中， `SequenceAlignment` 采用的也许是那种算法)。借助于 Mathematica 的前端， 这个函数将会有更漂亮的输出结果。(这里插一句，Unix 下面的 [`diff`][3] 也是用 LCS 算法实现的)。

PS:

*   在文章的前一部分讨论了部分实现的原理， 如果读者只是希望在 Mathematica 中有一个很好的 `diff` 函数的话， 可以跳过那些内容，直接复制源代码就可以了。也可以直接在[这里][4]直接下载我写好的包（package）。
*   要理解源代码， 读者至少需要对 Mathematica 中的模式匹配，函数式编程都有初步的了解。我不会在此文中解释 `/@` 是什么算符，`__String` 又能代表什么。读者可以查 Mathematica 文档去了解他们的用途.

`[wldoc]SequenceAlignment[/wldoc]` 是我们实现 `diff` 所需要用到的核心内置函数。出于效率的考虑， 在 Mathematica 的实现中， 它是直接用 `C` 写的。在我们的例子， 并不用考虑该[文档][2]中的诸多选项参数，我们只需要用到这个函数最简单的功能。`SequenceAlignment` 可以作用在两个表（list）上， 例如:

[wlcode] In[x] := SequenceAlignment[{1, 2, 3, 4, 5, 6}, {1, 10, 11, 4, 5, 6}] [/wlcode]

[wlcode] Out[x] = {{1}, {{2, 3}, {10, 11}}, {4, 5, 6}} [/wlcode]

也可以作用在字符串上， 例如:

[wlcode] In[x] := SequenceAlignment["Proud of being a Chinese", "proud of being a Japanese"] Out[x] = {{"P", "p"}, "roud of being a ", {"Chi", "Japa"}, "nese"} [/wlcode]

先来解释一下后一个例子里的输出。它实际上就是把两个字符串公共的部分拎出来， 例如`"roud of being a "`和`"nese"`， 然后再把不同的部分放在列表里， 例如`{"P", "p"}`和`{"Chi", "Japa"}`。这里有个并不trival的地方需要说明一下，所谓的**公共部分**其实是有相当的歧义的， 对很多例子而言， 从不同的方式看，可以得到不同的公共部分， 例如:

[wlcode] In[x] := SequenceAlignment["234", "3426"]

Out[x] = {{"2", ""}, "34", {"", "26"}} [/wlcode]

可见默认地`"34"`被当做了公共部分。但是，我们其实也可以把`"2"`当做公共部分呀，然后得到相应的结果`{{"", "34"}, "2", {"34", "26"}}`。确实是可以的，`SequenceAlignment` 由于采用了 LCS 算法， 所以相当于采取的选择标准是**`公共序列的总长度最大化`**， 这一点尤其适合用来做文本的差异分析，因为最大化的公共子序列确定了两个被对比文本的固定结构。

对于 `SequenceAlignment` 作用在两个表上的那个例子，实际上就是把着眼的最小单位从单个字符转移到列表中的单个元素，其他的逻辑都是完全一样的。

粗略地解释了 `SequenceAlignment` 的行为后，我们开始逐渐地朝我们的目标前进。

## 原理简析

Unix 下的 `diff` 是一个 *Compare files line by line.* 的工具。所谓 *line by line* 当然不是傻傻地一行一行地比较，而是指在基于 LCS 算法确定的文本固定结构后再一行一行比较，所以你会发现即使你在旧文件的开头加了一行得到新文件，即使新旧文件所有的行都错位， **diff** 也能迅速察觉到，因为添加一行并不会阻碍 LCS 算法确定文本固定结构。有的读者会问， 为什么要compare line by line 而不把新旧文件整个而当做字符串，对比两个字符串的差异， 那样不是更精细么? 问题的答案是**效率**。事实上Mathematica 系统中包含的一个让内部开发人员使用，但是也可以在发行版中访问的函数 ``System`Dump`showStringDiffs[str1_String， str2_String]`` 就是这么干的， 它可以输出很漂亮的结果， 但是字符串一长， 就会导致 Mathematica 卡死。让我们来做一个简单的分析:

考虑我们有 A、B 两个文件， 文件大概都在 n = 100000 个字符左右，着眼到单个字符， 因为 LCS的复杂度为 $O(n^2)$，所以可以想象到需要运行的时间。现在假设 A B 文件大约每 k = 20 个字符构成一行， 着眼到单个行进行比较， 复杂度变成了 $O((\frac{n}{k})^2)$， 这个提升是非常明显的。但是这样做不够精细， 它只会告诉你哪些行被改了(修改前后的行以**行对**的形式给出，一个行对中包含修改前后的两个行)，而行内哪些字符被改动了却不得而知， 一个很自然却非常关键的想法是**分别对每个行对再作一次 `diff`， 并且这次的 `diff` 左眼于单个字符**，这一步的复杂度是 $O(\frac{n}{k}k^2)$， 也即 $O(nk)$。这样一来，连被修改的行中哪些字符被修改了，我们也知道了.

*Compare files line by line* 有非常重要的意义，要知道 `diff` 用于很多现代[版本控制系统][5]当中，而其中经常要做的是对修改前后的源代码进行差异分析。对于源代码而言，通常有很好的结构， 不会像普通文本一样成百上千个字符才来一个换行符，而且源代码的修改前后通常也是保持原有的大体结构。所以 diff 往往可以高效的运用在上面.

## `DiffList`

讲了这么多，下面是我们要写的第一个函数 `DiffList`， 我们希望达到的效果是 `DiffList` 接收两 list 作为参数， 用 `SequenceAlignment` 作差异分析， 然后返回格式化后的分析结果.

下面是我写好的满足要求的代码:

[wlcode] Clear[DiffList]; DiffList[list1_List, list2_List] := Module[ {tmp}, Flatten[#, 1]&@ ( Switch[ #, {{}, \_List}, {"", #}& /@ Last[#], {\_List, {}}, {#, ""}& /@ First[#], {\_List, \_List}, tmp = Max[Length/@#]; Transpose@(PadRight[#,tmp,""]& /@ #), _List, {#, #}& /@ # ]& /@ SequenceAlignment[list1,list2] ) ] [/wlcode]

测试一下运行结果:

[wlcode] In[x] := DiffList[{"Programming", "is", "so", "interesting"}, {"Programming", "in", "Mathematica", "is", "interesting"}] Out[x] = {{"Programming", "Programming"}, {"", "in"}, {"", "Mathematica"}, {"is", "is"}, {"so", ""}, {"interesting", "interesting"}}. [/wlcode]

这个输出的可读性比较差，我们用 [`Grid`][6] 函数把它排版一下显示给读者。

    Programming    Programming
                       in
                    Mathematica
        is             is
        so  
    interesting     interesting
    

把我们的话题限制于参数是两个 list of strings (String当然可以被替换成其他类型)。可以看见， `DiffList` 函数的返回结果是`行对`的列表， 每个`行对`之中的两个 `string` 被认为是在两个list上处于统一位置的 `string`。当然了， 两个 list 中，有时对应的一个位置， 只有其中一个 list 会有元素， 那么我们就用 `""` 当做是另外一个 list 对应位置的元素.

下面我们来稍微解释一下 `DiffList` 的实现，尽管我认为这对我们理解整个程序并没有太大意义。`SequenceAlignment[list1,list2]` 生成了类如`{{"Programming"}, {{}, {"in", "Mathematica"}}, {"is"}， {{"so"}, {}}, {"interesting"}}` 这样的表结构， 表中的元素则分成下面的两类:

*   `{__String}` 型， 代表列表里面的元素都是 Common Sequence。
*   `{_List, _List}` 型， 代表有差异的一些元素， `差异`可以是被增加，被删除， 或者被修改。

接着， 在 `DiffList` 中有个由 `Switch` 构造的函数被 `Map` 到了 `SequenceAlignment[list1,list2]` 的结果上， 也就是分别作用在每个元素上。`Switch` 的内部就是一系列的匹配， 下面依次说明:

*   `{{}, _List}, {"", #}& /@ Last[#]` 中 `{}` 表示左边的字符串们缺失， 所以我们用 `""` 填充；
*   `{_List, {}}, {#, ""}& /@ First[#],` 右边的字符串们缺失；
*   `{_List, _List}, tmp = Max[Length/@#]; Transpose@(PadRight[#,tmp,""]& /@ #)` 尽管两边都是非空 list， 但是很有可能里面包含的元素并不相同， 较短的用 `""` 填充使得两个 list 等长；
*   `_List, {#, #}& /@ #` 最后剩下的就是单纯的 list of strings 的情况了， 此时它们属于 Common Sequence， 为了统一格式， 复制一份， 变成 `{a_List, a_List}` 的形式；
*   最后再用 `Flatten[#, 1]&` 去除 `level 1` 的多余的 `list head`。

## `TextDiff`

现在利用 `DiffList` 函数我们已经可以找到哪些行被修改了，通过判断`行对`中是否包含 `""`， 可以推断出是以何种方式修改的: 添加， 删除，或者替换。接着我们按照此前分析的，对每个`行对`再进行一次差异分析，可以让我们得到更精细化的结果，精确到具体是哪个字符被修改了。下面的代码中利用定义好的 `DiffList` 以及内置的 `SequenceAlignment`，构建了一个满足我们最终需要的函数。只有一个 `AddStyle`函数是还没有定义的， 我们将在文章的后面给出它的源代码和说明。

[wlcode] Clear[TextDiff]; mcTextDiff[oldText_String,newText_String] := Module[ {tmp,cmpStrings}, tmp = DiffList[StringSplit[oldText,"\n"], StringSplit[newText,"\n"]]; cmpStrings[s1_,s2_] := Replace[SequenceAlignment[s1,s2],a_String:>{a,a},1]; tmp = Composition[ (Row /@ Transpose[#])& /@ #&, (AddStyle /@ #&) /@ #&, (cmpStrings[First[#], Last[#]]& /@ #)& ][tmp]; Grid[tmp, ItemSize->Fit, Alignment->Left, Frame->All, FrameStyle->LightGray] ] [/wlcode]

大体说明一下这个函数的定义:

*   `tmp = DiffList[StringSplit[oldText,"\n"], StringSplit[newText,"\n"]];`：`StringSplit[text, "\n"]` 用于把 `String` 类型的 `text` 换行符为标志分割开来， 反回一个分割后元素的列表。列表元素中不再包含换行符。之后利用 `DiffList`函数返回对齐的差异分析列表， 返回的结果可以表示出哪些行被修改了
*   `cmpStrings[s1_,s2_] := Replace[SequenceAlignment[s1,s2],a_String:>{a，a},1];` 定义一个临时的函数, 返回对齐的对字符串的差异分析结果.
*   `tmp = Composition[..][tmp]` 在没有每一步的输出结果的情况下，这段代码非常不容易解释， 读者最好在Mathematica亲自运行一下， 依次观察`Composition` 中每一个函数作用后的输出结果， 这里`Composition` 相当于将数学中函数复合的操作。
*   `Grid[tmp, ItemSize->Fit, Alignment->Left, Frame->All, FrameStyle->LightGray]` 最后一行虽然有很多参数， 但是只是为了让输出的结果好看一点， 直接改成 `Grid[tmp]` 也可以。`Grid` 函数可以让结构齐整的列表结构以表格的形式输出， 利于阅读。

## `AddStyle`

下面我们来补充说明 `AddStyle`。读者只要记住这段代码是用来对输入的型如 `{_String, _String}` 添加一些输出的格式， 譬如颜色， 字体， 外围方框等等，不理解这段代码对理解 `TextDiff` 的工作机制也没有什么影响。我给出源代码:

[wlcode] Clear[AddStyle]; Options[AddStyle] = { "Style" -> { LightGreen, LightRed, {LightPurple, LightBlue} } }; AddStyle[list_List, OptionsPattern[]] := Module[ { keptColor = First@OptionValue["Style"], removedColor = OptionValue["Style"][[2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2]], changedColors = Last@OptionValue["Style"], removedStyle = {Bold, FontVariations -> {"StrikeThrough" -> True}} }, Switch[ list, {a_, a_}, list, {"", *}, {Framed[#, Background -> removedColor]& @ Style[Last[list]， Sequence@@removedStyle], Framed[#, Background -> keptColor]& @ Style[Last[list], Bold]}, {*, ""}, {Framed[#， Background -> keptColor]& @ Style[First[list], Bold], Framed[#, Background -> removedColor]& @ Style[First[list], Sequence@@removedStyle]}, {\_, \_}, {Framed[#, Background -> First@changedColors]& @ Style[First[list], Bold], Framed[#, Background -> Last@changedColors]& @ Style[Last[list], Bold]} ] ] [/wlcode]

上面的三段代码 `DiffList`， `TextDiff`， `AddStyle` 就完成了我们预订的目标， 这里只有 `TextDiff` 是供外部用户调用的， 其他两个都属于辅助性质的函数。现在我们测试一下函数的效果。我把 `new_LCS.cpp` 和 `old_LCS.cpp` 两个文件放在 `~/` 目录下， 由于 `TextDiff` 的参数是字符串， 所以我们用 Import[file， “Text”] 来导入上面的两个文件。

[wlcode] TextDiff[Import["~/old_LCS.cpp", "Text"], Import["~/new_LCS.cpp", "Text"]] // Timing [/wlcode]

后置 Timing 是为了观察一下耗时。下面是输出结果:

![illustration for textDiff][7]

运行这个函数大概花费 0.03 sec，对大文件的测试结果也令人相当满意。`new_LCS.cpp` 和 `old_LCS.cpp` 的文件可以在[这里][8]找到。

 [1]: http://zh.wikipedia.org/wiki/%E6%9C%80%E9%95%BF%E5%85%AC%E5%85%B1%E5%AD%90%E5%BA%8F%E5%88%97
 [2]: http://reference.wolfram.com/mathematica/ref/SequenceAlignment.zh.html
 [3]: http://en.wikipedia.org/wiki/Diff
 [4]: https://github.com/MathCraft/MathCraftAddOn/blob/master/MathCraftAddOn/MathCraftAddOn/Text.m
 [5]: http://baike.baidu.com/view/183136.htm
 [6]: http://reference.wolfram.com/mathematica/ref/Grid.zh.html
 [7]: http://www.mathcraft.org/posts/wp-content/uploads/2015/07/diff.gif
 [8]: https://github.com/MathCraft/MathCraftAddOn/tree/master/MathCraftAddOn/MathCraftAddOn/TestData
