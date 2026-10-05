// =============================================================
// 第 1 课：一个 process + 一个 channel
// 运行方式：nextflow run lessons/01_hello.nf
//
// 你会看到 SAY_HELLO 被执行了 3 次 —— 通道里有几个元素，
// process 就自动运行几次，这就是 Nextflow 的隐式并行。
// =============================================================

nextflow.enable.dsl = 2

// process 是最小执行单元：声明输入、输出、以及要执行的脚本
process SAY_HELLO {
    input:
    val name            // val: 接收一个普通值（字符串、数字等）

    output:
    stdout              // 把脚本的标准输出捕获为输出

    script:
    """
    echo "你好, ${name}! 欢迎学习 Nextflow"
    """
}

workflow {
    // Channel.of() 创建一个值通道，元素依次流入下游 process
    ch_names = Channel.of('小单', '小细', '小胞')

    // 调用 process 并把输出通道用 view() 打印出来
    SAY_HELLO(ch_names).view()
}
